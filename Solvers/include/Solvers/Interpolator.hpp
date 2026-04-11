// Interpolator.hpp
//
// 
//
// JAF 1/19/2026 

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H 

#include<memory>
#include<vector>
#include<utility> // std::pair
#include<Eigen/Core>

#include "ImplicitSolver.hpp"

namespace Solvers{ 

template<typename LHS_EXPR, 
          typename RHS_EXPR, 
          typename OSTEP_TUP, 
          template<typename L,typename R, typename O> SOLVER_T = Solvers::ImplicitSolver, 
          typename M = LinOps::SharedConstMesh1D
>
class Interpolator
{
  protected:
    // Type Defs ----------------------------- 
    struct VecSaveWrite
    {
      std::vector<Eigen::VectorXd>& m_vec; 
      VecSaveWrite()=delete; 
      VecSaveWrite(std::vector<Eigen::VectorXd>& v_init) : m_vec(v_init){}; 
      VecSaveWrite(const VecSaveWrite& other) : m_vec(other.m_vec){}; 
      void SaveSolution(Eigen::VectorXd&& sol) const { m_vec.emplace_back(sol); }; 
      void ConsumeLastSolution(Eigen::VectorXd&& sol) const { m_vec.emplace_back(sol); }; 
    }; 

    // Member Data ------------------------------
    std::vector<Eigen::VectorXd> m_data; 
    SOLVER_T<LHS_EXPR, RHS_EXPR, OSTEP_TUP> m_solver; 
    SolverArgs<M> m_args; 
    std::shared_ptr<const LinOps::MeshXD> m_high_dim_mesh; 
    bool m_calculated; 

  public: 
    // Constructors + Destructor ===================================
    // default 
    Interpolator()=delete; 

    // from args 
    Interpolator(LHS_EXPR& lhs, RHS_EXPR& rhs, OSTEP_TUP ostep_init, SolverArgs<M> args_init = {})
      : m_data(0), 
      m_args(std::move(args_init)), 
      m_high_dim_mesh(LinOps::make_MeshXD(m_args.domain_mesh_ptr)), 
      m_solver(lhs,rhs,ostep_init), 
      m_calculated(false) 
    {}

    // copy 
    Interpolator(const Interpolator& other)=delete; 

    // destructor 
    ~Interpolator()=default; 

    // Member Funcs ====================================================
    
    // get value of Solution at any point t,x1,x2,...xn in time/space 
    template<typename... Scalars>
    double SolAt(double t, Scalars... coords)
    {
      std::array<double, sizeof...(coords)> c{coords...};
      return SolAt(t, c); 
    }
    
    // specialization for 1D case 
    double SolAt(double t, double x)
    {
      std::array<double,1> c{x}; 
      return SolAt(t,c); 
    }
    
    // ... with container of spatial coords
    template<typename Cont_C>
    double SolAt(double t, const Cont_C& coords){
      // if m_data is empty... 
      if(!m_calculated) FillVals(); 

      // find index in m_args.time_mesh_ptr
      auto time_interval_pair = get_interval(*m_args.time_mesh_ptr, t);
      auto time_idx = std::distance(m_args.time_mesh_ptr->cbegin(), time_interval_pair.first); 

      // find left / right value in linear interpolation 
      double val_01 =  LinearInterp(coords, m_data[time_idx]); 
      double val_02 =  LinearInterp(coords, m_data[time_idx+1]);
      
      // linear interpolation (t-t1) * (y2 - y1) / (t2 - t1) 
      return val_01 + (t - *time_interval_pair.first) * (val_02 - val_01) / (*time_interval_pair.second - *time_interval_pair.first); 
    }

    // resets m_calculated to false. resize m_data
    void Reset()
    {      
      m_data.resize(0); 
      m_high_dim_mesh = LinOps::make_MeshXD(m_args.domain_mesh_ptr); 
      m_calculated=false; 
    }
    
    // getters to m_args
    const auto& Args() const { return m_args; }; 

    // getters to m_solver 
    auto& Solver(){return m_solver;} // non const allows access to solvers internals. i.e. SetMaxIterations 
    const auto& Solver() const {return m_solver;}
    
    // set m_args to a new input
    void SetArgs(SolverArgs<M> args_switch)
    {
      m_args = std::move(args_switch); 
      Reset(); 
    }

    // Getters to m_data 
    const auto& StoredData() const { return m_data; }; 

    // Populate m_data with solutions at each step in time 
    void FillVals()
    {
      if(!m_calculated)
      {
        // resize + reserve data 
        m_data.resize(0); 
        m_data.reserve(m_args.time_mesh_ptr->size());

        // WritePolicy moves all solutions at each time step to m_data
        m_solver.Calculate(m_args, VecSaveWrite(m_data)); 

        // update status of interp 
        m_calculated = true; 
      }
    }

  private:
    // Unreachable ----------------------------------------------------
    template<typename Cont_C, typename Cont_V>
    double LinearInterp(const Cont_C& coords, const Cont_V& v)
    {
      return LinearInterp_recursive_impl(coords, v, m_high_dim_mesh, coords.size()-1);
    }; 

    template<typename Cont_C, typename Cont_V>
    double LinearInterp_recursive_impl(
      const Cont_C& coords, 
      const Cont_V& v, 
      const std::shared_ptr<const LinOps::MeshXD>& m,
      std::size_t ith_dim,
      std::size_t cumulative_offset = 0
    )
    {
      const auto& sub_dim_m = m->getMesh1DSafe(ith_dim); 
      auto bounding_interval = get_interval(*sub_dim_m, coords[ith_dim]);  

      if(ith_dim == 0)
      {
        std::size_t final_offset_01 = cumulative_offset + std::distance(sub_dim_m->cbegin(), bounding_interval.first); 
        std::size_t final_offset_02 = final_offset_01 + 1; 
        // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
        double result = v[final_offset_01] + (coords[ith_dim] - *bounding_interval.first) * (v[final_offset_02]-v[final_offset_01]) / (*bounding_interval.second - *bounding_interval.first);  
        return result; 
      }
      else
      {
        std::size_t stride_size = m->sizesMiddleProduct(0,ith_dim); 
        std::size_t interval_start_idx = std::distance(sub_dim_m->cbegin(), bounding_interval.first); 
        std::size_t offset_01 = cumulative_offset + stride_size * (interval_start_idx); 
        std::size_t offset_02 = cumulative_offset + stride_size * (interval_start_idx+1); 

        double endpoint_01 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, offset_01); 
        double endpoint_02 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, offset_02);
        
        // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
        double result = endpoint_01 + (coords[ith_dim] - *bounding_interval.first) * (endpoint_02 - endpoint_01) / (*bounding_interval.second - *bounding_interval.first); 
        return result; 
      }
    }; 

    template<typename Cont>
    auto get_interval(const Cont& v, double c)
    {
    // runtime checks 
    if(v.size() < 2) throw std::runtime_error("size of v < 2"); 
    if(c < v.cbegin()[0]) throw std::runtime_error("c < v[0]"); 

    // right side a(i+1) 
    auto after = std::lower_bound(v.cbegin(), v.cend(), c);
    if(after == v.cend()) throw std::runtime_error("right bound == v.cend()"); 

    // if a(i+1) == v[0] bump it by 1. 
    auto before = (after==v.cbegin()) ? after++ : std::prev(after); 
    return std::pair(before, after); 
    }; 

}; 

} // end namespace Solvers 

#endif // Interpolator.hpp


// Linear interpolation scratch work 


// 1D Case: 

/* 
Find pair [x1,x2] with x1 <= x <= x2] and return (x-x1) * (val[idx(x2)] - val[idx(x1)]) / (x2-x1)
*/

// 2D Case:

/* 
find pair [y1,y2] with y1 <= y <= y2.

get idx(y1) and idx(y2)

along OneDim_view( idx(y1) ) 
get [x1,x2] with x1 <= x <= x2 and return ... 

along OneDim_view( idx(y2) )
get [x1,x2] ......

return (y-y1) * (1D_Interp(idx(y2)) - 1D_Interp(idx(y1))) / (y2-y1)
*/

// 3D Case:

/*
find pair [z1,z2] with z1 <= z <= z2 

get idx(z1) and idx(z2) 

along TwoDim_view( idx(z1) ) 
get [y1,y2] -> get [x1,x2] (twice) ... return 

along TwoDim_view( idx(z2) ) 
get [y1,y2] -> get [x1,x2] (twice) ... return 

return (z-z1) * (2D_Interp(idx(z2)) - 2D_INterp(idx(z1))) / (z2 - z1)

*/

// In General ND Case:

/*
find pair [x1, x2] with x1 <= x <= x2 


( 
  necessary data for N-1 Interpolation: 
    std::size_t N-1 : next inner dimension to get new [x1,x2] pairs in
    std::size_t k : offset to get solution value at in higher dimension   
)
along N-1 Dim View ( idx(x1) ) get val 1 
along N-1 Dim View ( idx(x2) ) get val 2

return (x - x1) * (val2 - val1) / (x2 - x1)
*/


/*
find pair [x1, x2] with x1 <= x <= x2 
auto p = get_interval(*m->getMesh1DSafe( current_ith_dimension ) )

if this is first dimension 
{
  ! requires offsets of all higher dimension. 
  ! actually return disc(idx(...) + some_offset_from_higher_dims)
  return discretization(idx(p.first)), discretization(idx(p.second)) 
}

if this is an intermediate dimension
{
  auto pair1 = LinearInterp_recursive_impl(...) with idx(p.first), (current dimension-1)
  auto pair2 = LinearInterp_recursive_impl(...) with idx(p.second), (current dimension-1)

  return (middle of pair1), (middle of pair2) 
}

( 
  necessary data for N-1 Interpolation: 
    std::size_t N-1 : next inner dimension to get new [x1,x2] pairs in
    std::size_t k : offset to get solution value at in higher dimension   
)
*/



// else if(ith_dim == (coords.size()-1)) // this is the outtermost dimension. 
// {
// // calculate new offset by incorporating this 
// std::size_t cumulative_offset_01 = cumulative_offset + m->sizesMiddleProduct(0, ith_dim) * std::distance(sub_dim_m->cbegin(), internval_pair.first); 
// std::size_t cumulative_offset_02 = cumulative_offset + m->sizesMiddleProduct(0, ith_dim) * std::distance(sub_dim_m->cbegin(), internval_pair.second); 

// // get 2 pairs of linear interpolant along subdimension (ith_dim-1)
// std::pair<double,double> val_pair_01 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, cumulative_offset_01);
// std::pair<double,double> val_pair_02 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, cumulative_offset_02);

// // perform linear interpolation in this dimension to produce pair of values 
// double val_01 = val_pair_01.first + (coords[ith_dim] - *internval_pair.first) * (val_pair_01.second-val_pair_01.first) / (*internval_pair.second-*internval_pair.first);
// double val_02 = val_pair_02.first + (coords[ith_dim] - *internval_pair.first) * (val_pair_02.second-val_pair_02.first) / (*internval_pair.second-*internval_pair.first);

// // return std::pair<double,double> 
// return {val_01, val_02}; 
// }