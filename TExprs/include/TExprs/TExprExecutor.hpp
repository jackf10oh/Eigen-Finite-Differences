// LhsExecutor.hpp
//
//
//
// JAF 1/16/2025 

#ifndef LHSEXECUTOR_H
#define LHSEXECUTOR_H 

#include<Eigen/Core>
// #include<LinOps/Mesh.hpp> // MeshPtr_t
// #include<LinOpsXD/MeshXD.hpp> // MeshPtr_t
#include<Utilities/SparseDiagExpr.hpp>
#include<Utilities/FornbergCalc.hpp> 
#include "TimeDerivBase.hpp" // MatrixStorage_t

namespace TExprs{

namespace traits{
// =============================================================================
template<typename TIMEDERIV_T>
struct TimeDerivTraits
{
  using _coeffatreturntype_unclean_ = decltype(std::declval<TIMEDERIV_T>().CoeffAt(std::declval<std::vector<double>>(),0,0));  
  using CoeffAtReturnType = std::remove_cv_t<std::remove_reference_t<_coeffatreturntype_unclean_>>;
}; 

template<typename TIMEDERIV_T>
struct coeffat_returns_double : public std::is_same<double, typename TimeDerivTraits<TIMEDERIV_T>::CoeffAtReturnType>{}; 

template<typename TIMEDERIV_T>
struct coeffat_returns_other : public std::negation<coeffat_returns_double<TIMEDERIV_T>>{}; 

} // end namespace traits 

// ===========================================================================
template<template<typename...> class PRED, typename TUP_T>
auto filter_tup(TUP_T tup)
{
  auto filter_lam = [](auto&& elem){
    using CLEAN_T = std::remove_reference_t<std::remove_cv_t<decltype(elem)>>;
    if constexpr(PRED<CLEAN_T>::value){
      // singleton tuple ( elem )
      return std::make_tuple(std::forward<decltype(elem)>(elem)); 
    }
    else{
      // empty tuple ( _ )
      return std::tuple<>{}; 
    }
  }; 

  auto singletons = std::apply(
    [&](auto&&... elems){
      return std::make_tuple(filter_lam(std::forward<decltype(elems)>(elems))...); 
    }, 
    tup
  ); 

  return std::apply(
    [](auto&&... elems){
      return std::tuple_cat(std::forward<decltype(elems)>(elems)...); 
    }, 
    singletons
  ); 
}

// ===============================================================================
template<typename TEXPR_T>
class TExprExecutor
{
  public:
    // Type Defs ------------------------------------------- 
    using TUP_T = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<TEXPR_T>().toTuple())>>;
    using SCALAR_TUP_T = std::remove_cv_t<std::remove_reference_t<decltype(filter_tup<TExprs::traits::coeffat_returns_double>(std::declval<TUP_T&>()))>>;
    using MAT_TUP_T = std::remove_cv_t<std::remove_reference_t<decltype(filter_tup<TExprs::traits::coeffat_returns_other>(std::declval<TUP_T&>()))>>;
    using INV_COEFF_T = std::conditional_t<std::tuple_size<MAT_TUP_T>::value==0, double, MatrixStorage_t>; 
  
  private:
    // Member Data ---------------------------------
    // tuple that stores all entries in expr_init.toTuple() such that .coeffAt(...) returns a double  
    SCALAR_TUP_T m_scalar_coeff_sum_partition; 
    // .....  such that .coeffAt(...) returns a Matrix  
    MAT_TUP_T m_mat_coeff_sum_partition; 

    // Highest order weights will be calculated to
    std::size_t m_order; 
    // number of nodes used in Fornberg algorithm 
    std::size_t m_num_nodes; 

    // list of t0, t1, ..., tn 
    std::vector<double> m_stored_times; 
    // list of solutions u0, u1, ..., un-1 at times t0, t1, ..., tn-1 
    std::vector<Eigen::VectorXd> m_stored_sols; 

    // forberg weights calculator 
    FornCalc m_weights_calc; 

    // result of BuildRhs 
    Eigen::VectorXd m_rhs_vec; 
    
    // result of build_inv_coeff 
    INV_COEFF_T m_inv_coeff; 

  public:
    // Constructors + Destructor =======================================
    // default 
    TExprExecutor()=delete; 
    // from expressions of Time derivatives
    TExprExecutor(TEXPR_T& expr_init) 
      :m_order(expr_init.Order()), 
      m_num_nodes(expr_init.Order()+1),  
      m_weights_calc(m_num_nodes, m_order), 
      m_stored_times(m_num_nodes), 
      m_stored_sols(m_order),
      m_scalar_coeff_sum_partition(filter_tup<TExprs::traits::coeffat_returns_double>(expr_init.toTuple())), 
      m_mat_coeff_sum_partition(filter_tup<TExprs::traits::coeffat_returns_other>(expr_init.toTuple()))
    {}
    // copy 
    TExprExecutor(const TExprExecutor& other)=delete; 
    // destructor 
    ~TExprExecutor()=default; 

    // Member Funcs =======================================
    // returns ref to newest solution 
    const auto& MostRecentSol() const { return m_stored_sols[m_stored_sols.size()-1]; }
    auto& MostRecentSol(){ return m_stored_sols[m_stored_sols.size()-1]; }

    // return ref to first elem in m_stored_sols. Gives an opportunity to move it elsewhere before overwritten in ConsumeSolution  
    Eigen::VectorXd& ExpiringSol(){ return m_stored_sols[0]; }

    // return full vector of StoredSols. 
    auto& StoredSols(){ return m_stored_sols; }
    
    // consume a time. push back all previous
    void ConsumeTime(double t)
    {
      // iterate from m_stored_times[0] ... [n-2]
      std::move(std::next(m_stored_times.begin()), m_stored_times.end(), m_stored_times.begin()); 
      m_stored_times[m_num_nodes-2] = t;
      // all values in m_stored_time have been left shifted by 1
      // first value dropped.
      // second from right most value == t 
      // right most value unassigned 
    }

    // ConsumeTime but copies full list into m_stored_times
    template<typename INPUT_IT>
    void ConsumeTimeList(INPUT_IT start, INPUT_IT end)
    {
      if(std::distance(start,end)!=m_num_nodes-1) throw std::runtime_error("Error: distance(start,end) must == size of m_stored_times - 1"); 
      std::move(
        std::move_iterator(start),
        std::move_iterator(end),
        m_stored_times.begin() 
      );
    }

    // consume a solution. push back all previous 
    void ConsumeSolution(Eigen::VectorXd sol)
    { 
      /* same idea as ConsumeTime(). with move semantics*/
      std::move(std::next(m_stored_sols.begin()), m_stored_sols.end(), m_stored_sols.begin()); 
      m_stored_sols[m_stored_sols.size()-1] = std::move(sol); 
    }

    // ConsumeSolution but copies full list into m_stored_sols
    template<typename INPUT_IT>
    void ConsumeSolutionList(INPUT_IT start, INPUT_IT end)
    {
      if(std::distance(start,end)!=m_stored_sols.size()) throw std::runtime_error("Error: distance(start,end) must == size of m_stored_sols"); 
      std::move(
        std::move_iterator(start),
        std::move_iterator(end),
        m_stored_sols.begin() 
      );
    }


    template<typename ANYMESHPTR_T>
    void set_mesh(ANYMESHPTR_T m)
    {
      auto set_mesh_lam = [&](auto&& elem){
        elem.set_mesh(m); 
      };
      auto set_mesh_for = [&](auto&&... elems){(set_mesh_lam(std::forward<decltype(elems)>(elems)),...);}; 
      // std::apply(set_mesh_for, m_scalar_coeff_sum_partition); // unnecessary. 
      std::apply(set_mesh_for, m_mat_coeff_sum_partition); 
    }

    void expr_SetTime(double t)
    {
      auto SetTime_lam = [&](auto&& elem){
        elem.setTime(t); 
      };
      auto set_mesh_for = [&](auto&&... elems){(SetTime_lam(std::forward<decltype(elems)>(elems)),...);}; 
      // std::apply(set_mesh_for, m_scalar_coeff_sum_partition); // unnecessary.  
      std::apply(set_mesh_for, m_mat_coeff_sum_partition); 
    }
    
    // Builds m_rhs_vector + m_inv_coeff from the next time; 
    void BuildNextTime(double next_t)
    {
      // update m_weights_calc to new time step.
      SetWeightsFromTime(next_t); 
      // ... assemble m_inv_coeff 
      m_inv_coeff = build_inv_coeff();
      // use m_weights_calc to assemble RHS
      m_rhs_vec = BuildRhs(); 
    }

    // getters to m_inv_coeff + m_rhs_vec 
    const INV_COEFF_T& inv_coeff() const { return m_inv_coeff; }
    Eigen::VectorXd& RhsVector(){ return m_rhs_vec; } 

  private:
    // Unreachable =========================================== 
    // from a time. set weights
    void SetWeightsFromTime(double t)
    {
      // set last entry in m_stored_times
      m_stored_times[m_num_nodes-1] = t; 
      // recalculate m_weights_calc 
      m_weights_calc.Calculate(t, m_stored_times.cbegin(), m_stored_times.cend(), m_order);     
    } // end SetWeightsFromTime 

    // // Build RHS starting point from m_stored_sols[0, ..., N-2]
    Eigen::VectorXd BuildRhs()
    {
      // initialize RHS vector to all zeros 
      Eigen::VectorXd result(m_stored_sols[0].size());
      result.setZero(); 

      // starting at oldest solution, first node's weight 
      for(std::size_t ith_node = 0; ith_node<m_num_nodes-1; ith_node++)
      {
        // if there are any scalar coeffs start with those 
        if constexpr(std::tuple_size<SCALAR_TUP_T>::value > 0)
        {
          double s = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, ith_node) + ...); 
            }, 
            m_scalar_coeff_sum_partition
          ); 
          // std::cout << "BuildRhs: " << s << std::endl;
          result += s * m_stored_sols[ith_node]; 
        } 
        if constexpr(std::tuple_size<MAT_TUP_T>::value > 0)
        {
          auto M = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, ith_node) + ...); 
            }, 
            m_mat_coeff_sum_partition 
          ); 
          // std::cout << "BuildRhs Mat: " << M << std::endl;
          result += M * m_stored_sols[ith_node]; 
        }
      } // end for() over ith_node 

      // flip result by negative 1 
      result *= (-1.0); 

      // multiply by inv_coeff_util; 
      if constexpr(std::tuple_size<MAT_TUP_T>::value == 0){
        // INV_COEFF_T is a scalar
        result *= m_inv_coeff; 
      }
      else{
        // INV_COEFF_T is a Matrix 
        Eigen::VectorXd temp = m_inv_coeff * result; 
        result = std::move(temp); 
      }
      // Eigen::VectorXd 
      return result; 
    }

    // gets 1 / c where c is coeff of U(n+1) in fdm equation 
    auto build_inv_coeff()
    {
      // all CoeffAt's evaluate to scalar -> return 1 / sum(coeffs...)
      if constexpr(std::tuple_size<MAT_TUP_T>::value == 0){
        double s = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, m_num_nodes-1) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        return 1.0 / s;  
      }
      else if constexpr(std::tuple_size<SCALAR_TUP_T>::value == 0){ 
        // all COeffAt's evaluate to matrix -> return (sum(coeffs)).cwiseInverse 
        MatrixStorage_t A = std::apply(
              [&](auto&&... coeffs){
                return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, m_num_nodes-1) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 
        MatrixStorage_t B = A.cwiseInverse(); 
        return B; 
      }
      else{ 
        // throw std::runtime_error("This formula hasn't been implemented yet!"); 
        // otherwise return product of 1/sum(scalar) + (sum(Mats)).cwiseInverse()
        double s = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, m_num_nodes-1) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        MatrixStorage_t A = std::apply(
              [&](auto&&... coeffs){
                return (coeffs.CoeffAt(m_weights_calc.m_arr, m_num_nodes, m_num_nodes-1) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 

        // combine s, A along diagonal 
        A = A.cwiseInverse(); 
        s = (1.0) / s; 
        for(std::size_t i=0; i<A.rows(); ++i) A.coeffRef(i,i) += s; 
        return A; 
      }
    } // end inv_coeff_util() 

}; 

} // end namespace TExprs  

#endif // LhsExecutor.hpp