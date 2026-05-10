// Interpolator.hpp
//
// Lazy calculator that uses a one of fdm's solvers 
// to write solutions at each time step into a stored vector
// and later interpolate any point  in time + space 
//
// JAF 4/12/2026 

#ifndef SOLUTIONINTERPOLATOR_H
#define SOLUTIONINTERPOLATOR_H 

#include<cassert>
#include<memory>
#include<vector>
#include<Eigen/Core>
#include "../Mesh.hpp"
#include "../LinOps/LinOpTraits.hpp" 
#include "../Utilities/LinearInterpolation.hpp"

#include "ImplicitSolver.hpp"

namespace fdm{
  namespace solvers{ 

template<
  typename StoredSolver, 
  typename C = Eigen::VectorXd // container of times 
>
class Interpolator
{
  protected:
    // Type Defs ----------------------------- 
    struct BackInserterSaver
    {
      std::vector<fdm::Vector>& m_vec; 
      BackInserterSaver()=delete; 
      BackInserterSaver(std::vector<Eigen::VectorXd>& v_init) : m_vec(v_init){}; 
      void saveSolution(Eigen::VectorXd sol){ m_vec.emplace_back(std::move(sol)); }; 
      void saveLastSolution(Eigen::VectorXd sol){ m_vec.emplace_back(std::move(sol)); }; 
    }; 

    // Member Data ------------------------------
    std::vector<Eigen::VectorXd> m_data; 
    StoredSolver m_solver; 
    SolverArgs<const fdm::Mesh,C> m_args; 
    bool m_calculated; 

  public: 
    // Constructors + Destructor ===================================

    // no default 
    Interpolator()=delete; 

    // from solver + args 
    Interpolator(StoredSolver s, SolverArgs<const fdm::Mesh,const C> args = {})
      : m_data(0), 
      m_args(std::move(args)), 
      m_solver(std::move(s)),
      m_calculated(false) 
    {}

    // not copyable! 
    Interpolator(const Interpolator& other)=delete; 

    // moveable 
    Interpolator(Interpolator&& other)=default; 

    // destructor 
    ~Interpolator()=default; 

    // Member Funcs ====================================================

    // getter to m_calculated
    inline bool isCalculated() const { return m_calculated; }

    // resets m_calculated to false. resize m_data
    void clearStoredSolutions()
    {      
      m_data.clear(); 
      m_calculated=false; 
    }
    
    // getters to m_args
    const auto& getArgs() const { return m_args; }; 

    // set m_args to a new input
    void setArgs(SolverArgs<const fdm::Mesh,const C> args_switch)
    {
      m_args = std::move(args_switch); 
      clearStoredSolutions(); 
    }

    // getters to m_solver 
    auto& getSolver(){return m_solver;} // non const allows access to solvers internals. i.e. SetMaxIterations 
    const auto& getSolver() const {return m_solver;}
    
    // Getters to m_data 
    const auto& getStoredSolutions() const { return m_data; }; 

    // Populate m_data with solutions at each step in time 
    void fillStoredSolutions()
    {
      if(m_calculated) return; 
      // resize + reserve data 
      m_data.resize(0); 
      m_data.reserve(m_args.times->size());

      // WritePolicy moves all solutions at each time step to m_data
      m_solver.calculate(m_args, BackInserterSaver(m_data)); 

      // update status of interp 
      m_calculated = true; 
    }

    // get value of Solution at any point t,{x1,x2,...xn} in time/space 
    template<std::size_t numDimsMax>
    double solAt(double t, const fdm::Coordinate<numDimsMax>& coords){
      assert((numDimsMax == m_args.mesh->numDims()) && "error in Interpolator.solAt() : # of dims in coordinate must == # of dims in stored mesh");
      // if m_data is empty... 
      if(!m_calculated) fillStoredSolutions(); 

      // find index in m_args.times
      auto time_interval = fdm::utils::make_subinterval(t, m_args.times->cbegin(), m_args.times->cend());
      auto offset = std::distance(m_args.times->cbegin(), time_interval.first); 

      // find left / right value in linear interpolation 
      double y1 =  interpolateSolution(coords.values, m_data[offset].cbegin(), m_data[offset].cend()); 
      double y2 =  interpolateSolution(coords.values, m_data[offset+1].cbegin(), m_data[offset+1].cend());
      
      // linear interpolation (t-t1) * (y2 - y1) / (t2 - t1) 
      return y1 + (t - *time_interval.first) * (y2 - y1) / (*time_interval.second - *time_interval.first); 
    }

  private:
    // Unreachable ----------------------------------------------------
    template<typename Container, typename Iterator>
    double interpolateSolution(const Container& coords, Iterator start, Iterator stop)
    {
      return LinearInterp_recursive_impl(coords, start, stop, m_args.mesh, coords.size()-1);
    }; 

    template<typename Container, typename Iterator>
    double LinearInterp_recursive_impl(
      const Container& coords, 
      Iterator start, 
      Iterator stop,  
      const std::shared_ptr<const fdm::Mesh>& m,
      std::size_t ith_dim,
      std::size_t cumulative_offset = 0)
    {
      const auto& sub_dim_m = m->getAxis(ith_dim); 
      auto subinterval = fdm::utils::make_subinterval(coords[ith_dim], sub_dim_m.cbegin(), sub_dim_m.cend());  

      if(ith_dim == 0){
        std::size_t final_offset = cumulative_offset + std::distance(sub_dim_m.cbegin(), subinterval.first); 
        auto it = std::next(start, final_offset); 
        double y1 = *it; ++it; 
        double y2 = *it; 
        // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
        return y1 + (y2-y1) * (coords[ith_dim] - *subinterval.first) / (*subinterval.second - *subinterval.first);  
      }
      else{
        std::size_t stride = m->sizesMiddleProduct(0,ith_dim); 
        std::size_t idx = std::distance(sub_dim_m.cbegin(), subinterval.first); 
        std::size_t next_offset = cumulative_offset + stride * (idx); 
        double y1 = LinearInterp_recursive_impl(coords, start, stop, m, ith_dim-1, next_offset); 
        double y2 = LinearInterp_recursive_impl(coords, start, stop, m, ith_dim-1, next_offset + stride);
        // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
        return y1 + (y2 - y1) * (coords[ith_dim] - *subinterval.first) / (*subinterval.second - *subinterval.first); 
      }
    }; 

}; 

  } // end namespace solvers
} // end namespace fdm 

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
// double y1 = val_pair_01.first + (coords[ith_dim] - *internval_pair.first) * (val_pair_01.second-val_pair_01.first) / (*internval_pair.second-*internval_pair.first);
// double y2 = val_pair_02.first + (coords[ith_dim] - *internval_pair.first) * (val_pair_02.second-val_pair_02.first) / (*internval_pair.second-*internval_pair.first);

// // return std::pair<double,double> 
// return {y1, y2}; 
// }