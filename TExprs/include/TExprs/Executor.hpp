// LhsExecutor.hpp
//
//
//
// JAF 1/16/2025 

#ifndef LHSTEXPREXECUTOR_H
#define LHSTEXPREXECUTOR_H 

#include<Eigen/Core>
#include<Utilities/SparseDiagExpr.hpp>
#include<Utilities/FornbergArrayCalc.hpp> 
#include "TExprTraits.hpp" 
#include "TimeDerivBase.hpp" // Matrix type

namespace TExprs{ 

// ===============================================================================
template<typename TimeDeriv, std::size_t M=TimeDeriv::maxOrder+1>
class Executor
{
  public:
    // Type Defs ------------------------------------------- 
    using Tup = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<TimeDeriv&>().toTuple())>>;
    using ScalarTup = std::remove_cv_t<std::remove_reference_t<decltype(TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_double>(std::declval<Tup>()))>>;
    using MatrixTup = std::remove_cv_t<std::remove_reference_t<decltype(TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_other>(std::declval<Tup>()))>>;
    using InvCoeff = std::conditional_t<std::tuple_size<MatrixTup>::value==0, double, TExprs::Matrix>; 

    // number of nodes used in Fornberg algorithm 
    static constexpr std::size_t numNodes = std::max(M, TimeDeriv::maxOrder+1); 

  private:
    // Member Data ---------------------------------

    // tuple that stores all entries in expr_init.toTuple() such that .coeffAt(...) returns a double  
    ScalarTup m_scalar_coeff_sum_partition;

    // .....  such that .coeffAt(...) returns a Matrix  
    MatrixTup m_mat_coeff_sum_partition; 

    // list of t0, t1, ..., tn 
    std::array<double, numNodes> m_stored_times; 

    // list of solutions u0, u1, ..., un-1 at times t0, t1, ..., tn-1 
    std::array<Eigen::VectorXd, numNodes-1> m_stored_sols; 

    // forberg weights calculator 
    FornArrayCalc<numNodes, TimeDeriv::maxOrder> m_weights_calc; 
    
    // result of build_inv_coeff.   
    InvCoeff m_inv_coeff; 

    // reference to the expression the executor is working on. 
    TimeDeriv& m_wrapped;

  public:
    // Constructors + Destructor =======================================
    // default 
    Executor()=delete; 

    // from expressions of Time derivatives
    Executor(TimeDeriv& expr_init) 
      : m_wrapped(expr_init), 
      m_scalar_coeff_sum_partition(TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_double>(m_wrapped.toTuple())), 
      m_mat_coeff_sum_partition(TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_other>(m_wrapped.toTuple()))
    {
      static_assert(TExprs::traits::is_timederiv_crtp<TimeDeriv>::value, "Must construct Executor from TExpr!"); 
    }

    // copy 
    Executor(const Executor& other)=default; 

    // move 
    // Executor(Executor&& other)=default; 
    // destructor 
    ~Executor()=default; 

    // Member Funcs =======================================

    // Getters to wrapped expression
    TimeDeriv& getWrappedExpression(){ return m_wrapped; } 
    const TimeDeriv& getWrappedExpression() const { return m_wrapped; } 

    // Getters to stored times. 
    auto& getStoredTimes(){ return m_stored_times; }; 
    const auto& getStoredTimes() const { return m_stored_times; }; 

    // getter to time that matches last solution in list 
    double& getCurrentTime(){ return m_stored_times[numNodes-2]; }
    const double& getCurrentTime() const { return m_stored_times[numNodes-2]; }

    // getter to future time
    double& getNextTime(){ return m_stored_times[numNodes-1]; }
    const double& getNextTime() const { return m_stored_times[numNodes-1]; }

    // getters to stored solutions 
    auto& getStoredSolutions(){ return m_stored_sols; }
    const auto& getStoredSolutions() const { return m_stored_sols; }

    // returns ref to newest solution 
    Eigen::VectorXd& getCurrentSolution(){ return m_stored_sols.back(); }
    const Eigen::VectorXd& getCurrentSolution() const { return m_stored_sols.back(); }

    // return ref to first elem in m_stored_sols. Gives an opportunity to move it elsewhere before overwritten in ConsumeSolution  
    Eigen::VectorXd& getExpiringSolution(){ return m_stored_sols.front(); }
    const Eigen::VectorXd& getExpiringSolution() const { return m_stored_sols.front(); }
    
    // consume a time. push back all previous
    void pushTime(double t)
    {
      // all values in m_stored_time have been left shifted by 1
      std::move(std::next(m_stored_times.begin()), m_stored_times.end(), m_stored_times.begin()); 
      // most recent time set to t 
      m_stored_times.back() = t;
    }

    // sets current times from start,end iterators.  
    template<typename Iter>
    void pushTimeRange(Iter start, Iter end){ 
      // if std::distance(start,end) != numNodes, 
      // only copy to rightmost times 
      std::size_t d = std::distance(start,end); 
      if(d > numNodes) throw std::runtime_error("Executor setStoredTimes(it,it) error: distance(start,end) > numNodes"); 
      std::copy(std::prev(m_stored_times.end(),numNodes-d), m_stored_times.end(), m_stored_times.begin()); 
      std::copy(start,end,std::next(m_stored_times.begin(),numNodes-d)); 
    }

    // consume a solution. push back all previous 
    void pushSolution(Eigen::VectorXd sol)
    { 
      /* same idea as ConsumeTime(). with move semantics*/
      std::move(std::next(m_stored_sols.begin()), m_stored_sols.end(), m_stored_sols.begin()); 
      m_stored_sols.back() = std::move(sol); 
    }

    // ConsumeSolution but copies full list into m_stored_sols
    template<typename Iter>
    void pushSolutionRange(Iter start, Iter end)
    {
      std::size_t d = std::distance(start,end); 
      if(d > numNodes-1) throw std::runtime_error("Executor pushSolutionRange(it,it) error: distance(start,end) > numNodes"); 
      // push numNodes - 1 - d solutions to front of m_stored_sols 
      std::move(
          std::move_iterator(std::prev(m_stored_sols.end(), d)),
          std::move_iterator(m_stored_sols.end()),
          m_stored_sols.begin()
      ); 
      // move [start,end) to last entries of m_stored_sols
      std::move(
          std::move_iterator(start),
          std::move_iterator(end),
          std::next(m_stored_sols.begin(), numNodes - 1 - d)
      ); 
    }

    // traverse stored tuples and sets mesh for any LinOps inside 
    template<typename AnyMeshSPtr>
    void setMesh(const AnyMeshSPtr& m)
    {
      std::apply(
        [&](auto&... args){ ((setMesh_singleton(m, args)), ...); }, 
        m_scalar_coeff_sum_partition
      ); 
      std::apply(
        [&](auto&... args){ ((setMesh_singleton(m, args)), ...); }, 
        m_mat_coeff_sum_partition
      ); 
    }

    // traverse stored tuples and sets time for any LinOps inside 
    void setTime(double t)
    {
      std::apply(
        [&](auto&... args){ ((setTime_singleton(t, args)),...); }, 
        m_scalar_coeff_sum_partition
      ); 
      std::apply(
        [&](auto&... args){ ((setTime_singleton(t, args)),...); }, 
        m_mat_coeff_sum_partition
      ); 
    }

    void calculate(){ 
      m_weights_calc.calculate(); 
      m_inv_coeff = buildInvCoeff(); 
    }

    decltype(auto) getRhsExpression()
    {
      return getRhsExpression_impl(std::make_index_sequence<numNodes-1>{}); 
    }

    // getters to m_inv_coeff + m_rhs_vec 
    const InvCoeff& getInvCoeff() const { return m_inv_coeff; }

  private:
    // Unreachable ===========================================     
    template<std::size_t... idxs>
    auto getRhsExpression_impl(std::index_sequence<idxs...>)
    {
      // assuming all the time derivatives don't have a time dependent CoeffOp, 
      // we can just use a super fast fold expression. 
      return -m_inv_coeff * (getRhsExpression_impl_helper<idxs> + ...); 
    }

    template<std::size_t ithNode>
    auto getRhsExpression_impl_helper()
    {
      if constexpr(std::tuple_size<MatrixTup>::value == 0){
        // just scalars need to be added 
        double scalar_sum = std::apply(
          [&](auto&&... args){
            return (args.template coeffAt<ithNode, numNodes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_scalar_coeff_sum_partition
        ); 
        // need to flip the coefficient by (-1) so that its moved to rhs 
        return scalar_sum * m_stored_sols[ithNode]; 
      }
      else if constexpr(std::tuple_size<ScalarTup>::value == 0){
        // just matrix return types need to be added
        decltype(auto) matrix_sum = std::apply(
          [&](auto&&... args){
            return (args.template coeffAt<ithNode, numNodes>(m_weights_calc.getArray()) + ...); 
          },
          m_mat_coeff_sum_partition
        ); 
        return matrix_sum * m_stored_sols[ithNode]; 
      }
      else{
        // sum of both scalar AND matrix return types! 
        double scalar_sum = std::apply(
          [&](auto&&... args){
            return (args.template coeffAt<ithNode, numNodes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_scalar_coeff_sum_partition
        );
        decltype(auto) matrix_sum = std::apply(
          [&](auto&&... args){
            return (args.template coeffAt<ithNode, numNodes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_mat_coeff_sum_partition
        ); 
        return scalar_sum * m_stored_sols[ithNode] + matrix_sum * m_stored_sols[ithNode]; 
      }
    }

    // gets 1 / c where c is coeff of U(n+1) in fdm equation 
    auto buildInvCoeff()
    {
      // all CoeffAt's evaluate to scalar -> return 1 / sum(coeffs...)
      if constexpr(std::tuple_size<MatrixTup>::value == 0){
        double s = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.template coeffAt<numNodes-1, numNodes>(m_weights_calc.getArray()) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        return 1.0 / s;  
      }
      else if constexpr(std::tuple_size<ScalarTup>::value == 0){ 
        // all COeffAt's evaluate to matrix -> return (sum(coeffs)).cwiseInverse 
        TExprs::Matrix A = std::apply(
              [&](auto&&... coeffs){
                return (coeffs.template coeffAt<numNodes-1, numNodes>(m_weights_calc.getArray()) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 
        TExprs::Matrix B = A.cwiseInverse(); 
        return B; 
      }
      else{ 
        // throw std::runtime_error("This formula hasn't been implemented yet!"); 
        // otherwise return product of 1/sum(scalar) + (sum(Mats)).cwiseInverse()
        double s = std::apply(
            [&](auto&&... coeffs){
              return (coeffs.template coeffAt<numNodes-1, numNodes>(m_weights_calc.getArray()) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        TExprs::Matrix A = std::apply(
              [&](auto&&... coeffs){
                return (coeffs.template coeffAt<numNodes-1, numNodes>(m_weights_calc.getArray()) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 

        // combine s, A along diagonal 
        A = A.cwiseInverse(); 
        s = (1.0) / s; 
        for(std::size_t i=0; i<A.rows(); ++i) A.coeffRef(i,i) += s; 
        return A; 
      }
    } 

    // Sets mesh onto a Time Derivative. traverse lhs/rhs of multiply expressions
    template<typename AnyMeshSPtr,typename TimeDerivU>
    void setMesh_singleton(const AnyMeshSPtr& m,TimeDerivU& tderiv)
    {
      if constexpr(TExprs::traits::is_coeffmult_crtp<TimeDerivU>::value){
        setMesh_singleton(m, tderiv.getLhs()); 
        setMesh_singleton(m, tderiv.getRhs()); 
      }
      if constexpr(LinOps::traits::is_linop_crtp<TimeDerivU>::value){
        tderiv.set_mesh(m); 
      }
    }

    // Sets mesh onto a Time Derivative. traverse lhs/rhs of multiply expressions
    template<typename TimeDerivU>
    void setTime_singleton(double t,TimeDerivU& tderiv)
    {
      if constexpr(TExprs::traits::is_coeffmult_crtp<TimeDerivU>::value){
        setTime_singleton(t, tderiv.getLhs()); 
        setTime_singleton(t, tderiv.getRhs()); 
      }
      if constexpr(LinOps::traits::is_linop_crtp<TimeDerivU>::value){
        tderiv.SetTime(t); 
      }
    }
}; 

// entry point for giving specific # of nodes, CTAD time derivate
template<std::size_t numNodes, typename TimeDeriv>
auto make_Executor(TimeDeriv& tderiv)
{
  return Executor<TimeDeriv,numNodes>(tderiv); 
} 

// entry point  without nodes being specified 
template<typename TimeDeriv>
auto make_Executor(TimeDeriv& tderiv)
{
  return Executor<TimeDeriv, TimeDeriv::maxOrder+1>(tderiv); 
}

} // end namespace TExprs  

#endif // LhsExecutor.hpp