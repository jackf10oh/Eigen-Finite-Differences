// LhsExecutor.hpp
//
//
//
// JAF 1/16/2025 

#ifndef LHSTEXPREXECUTOR_H
#define LHSTEXPREXECUTOR_H 

#include<algorithm>
#include<array>
#include<Eigen/Core>
#include "../Utilities/FornbergArrayCalc.hpp"
#include "TExprTraits.hpp" 
#include "../Types.hpp" // Scalar, DiagMatrix

namespace fdm{
  namespace texprs{

// ===============================================================================
template<typename TimeDeriv, std::size_t M=TimeDeriv::maxOrder+1>
class Executor
{
  public:
    // Type Defs ------------------------------------------- 
    using Tup = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<TimeDeriv&>().toTuple())>>;
    using ScalarTup = std::remove_cv_t<std::remove_reference_t<decltype(texprs::traits::filter_tup<texprs::traits::coeffat_returns_double>(std::declval<Tup>()))>>;
    using MatrixTup = std::remove_cv_t<std::remove_reference_t<decltype(texprs::traits::filter_tup<texprs::traits::coeffat_returns_other>(std::declval<Tup>()))>>;
    using InvCoeff = std::conditional_t<std::tuple_size<MatrixTup>::value==0, fdm::Scalar, fdm::DiagMatrix>; 

    // number of nodes used in Fornberg algorithm 
    static constexpr std::size_t numStoredTimes = std::max(M, TimeDeriv::maxOrder+1); 
    static constexpr std::size_t numStoredSols = numStoredTimes-1; 

  private:
    // Member Data ---------------------------------
    // reference to the expression the executor is working on. 
    TimeDeriv& m_wrapped;

    // tuple that stores all entries in expr_init.toTuple() such that .coeffAt(...) returns a double  
    ScalarTup m_scalar_coeff_sum_partition;

    // .....  such that .coeffAt(...) returns a Matrix  
    MatrixTup m_mat_coeff_sum_partition; 

    // list of t0, t1, ..., tn 
    std::array<double, numStoredTimes> m_stored_times; 

    // list of solutions u0, u1, ..., un-1 at times t0, t1, ..., tn-1 
    std::array<fdm::Vector, numStoredTimes-1> m_stored_sols; 

    // forberg weights calculator 
    fdm::utils::FornArrayCalc<numStoredTimes, TimeDeriv::maxOrder> m_weights_calc; 
    
    // result of buildInvCoeff.   
    InvCoeff m_inv_coeff; 

  public:
    // Constructors + Destructor =======================================
    // default 
    Executor()=delete; 

    // from expressions of Time derivatives
    Executor(TimeDeriv& expr_init) 
      : m_wrapped(expr_init),
      m_scalar_coeff_sum_partition(texprs::traits::filter_tup<texprs::traits::coeffat_returns_double>(m_wrapped.toTuple())), 
      m_mat_coeff_sum_partition(texprs::traits::filter_tup<texprs::traits::coeffat_returns_other>(m_wrapped.toTuple()))
    {
      static_assert(texprs::traits::is_timederiv_crtp<TimeDeriv>::value, "Must construct Executor from TExpr!");
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
    double& getCurrentTime(){ return m_stored_times[numStoredTimes-2]; }
    const double& getCurrentTime() const { return m_stored_times[numStoredTimes-2]; }

    // getter to future time
    double& getNextTime(){ return m_stored_times[numStoredTimes-1]; }
    const double& getNextTime() const { return m_stored_times[numStoredTimes-1]; }

    // getters to stored solutions 
    auto& getStoredSolutions(){ return m_stored_sols; }
    const auto& getStoredSolutions() const { return m_stored_sols; }

    // returns ref to newest solution 
    fdm::Vector& getCurrentSolution(){ return m_stored_sols.back(); }
    const fdm::Vector& getCurrentSolution() const { return m_stored_sols.back(); }

    // return ref to first elem in m_stored_sols. Gives an opportunity to move it elsewhere before overwritten in ConsumeSolution  
    fdm::Vector& getExpiringSolution(){ return m_stored_sols.front(); }
    const fdm::Vector& getExpiringSolution() const { return m_stored_sols.front(); }
    
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
      std::size_t d = std::distance(start,end);        
      if(d >= numStoredTimes) 
      {
        // throw std::runtime_error("Executor setStoredTimes(it,it) error: distance(start,end) > numStoredTimes");
        // only copy from rightmost times 
        // std::copy(std::next(m_stored_times.begin(),d-numStoredTimes), m_stored_times.end(), m_stored_times.begin()); 
        std::copy(std::next(start, d - numStoredTimes),end,m_stored_times.begin()); 
      } 
      else
      {
        // only copy into rightmost times 
        std::copy(std::next(m_stored_times.begin(),d), m_stored_times.end(), m_stored_times.begin()); 
        std::copy(start,end,std::next(m_stored_times.begin(),numStoredTimes-d)); 
      }
    }

    // rotates stored solutions. meaning the put to back of array. useful to keep memory allocated in hot loops
    void rotateStoredSolutions(std::size_t idx=1)
    {
      std::rotate(m_stored_sols.begin(), std::next(m_stored_sols.begin(),idx), m_stored_sols.end()); 
    }
    // consume a solution. push back all previous 
    void pushSolution(fdm::Vector sol)
    { 
      /* same idea as pushTime(). with move semantics*/
      std::move(std::next(m_stored_sols.begin()), m_stored_sols.end(), m_stored_sols.begin()); 
      m_stored_sols[numStoredSols-1] = std::move(sol); 
    }

    // ConsumeSolution but copies full list into m_stored_sols
    template<typename Iter>
    void pushSolutionRange(Iter start, Iter end)
    {
      std::size_t d = std::distance(start,end); 
      if(d >= numStoredTimes-1) 
      {
        // throw std::runtime_error("Executor pushSolutionRange(it,it) error: distance(start,end) > numStoredTimes"); 
        // move [start+d-numStoredTimes-1,end) to last entries of m_stored_sols
        std::move(
            std::move_iterator(std::next(start, d - (numStoredTimes-1))),
            std::move_iterator(end),
            m_stored_sols.begin()
        ); 
      }
      else
      {
        // push numStoredTimes - 1 - d solutions to front of m_stored_sols
        std::move(
            std::move_iterator(std::prev(m_stored_sols.end(), d)),
            std::move_iterator(m_stored_sols.end()),
            m_stored_sols.begin()
        ); 
        // move [start,end) to last entries of m_stored_sols
        std::move(
            std::move_iterator(start),
            std::move_iterator(end),
            std::next(m_stored_sols.begin(), numStoredTimes - 1 - d)
        ); 
      }
    }

    // traverse stored tuples and sets mesh for any LinOps inside 
    void setMesh(const std::shared_ptr<const fdm::Mesh>& m)
    {
      std::apply(
        [&m, this](auto&... args){ ((setMesh_singleton(m, args)), ...); }, 
        m_scalar_coeff_sum_partition
      ); 
      std::apply(
        [&m, this](auto&... args){ ((setMesh_singleton(m, args)), ...); }, 
        m_mat_coeff_sum_partition
      ); 
    }

    // traverse stored tuples and sets time for any LinOps inside 
    void setTime(double t)
    {
      std::apply(
        [t,this](auto&... args){ ((setTime_singleton(t, args)),...); }, 
        m_scalar_coeff_sum_partition
      ); 
      std::apply(
        [t,this](auto&... args){ ((setTime_singleton(t, args)),...); }, 
        m_mat_coeff_sum_partition
      ); 
    }

    void calculate(double t){ 
      m_weights_calc.calculate(t, m_stored_times.cbegin(), m_stored_times.cend()); 
      setTime(t);   
      buildInvCoeff(); 
    }

    decltype(auto) getRhsExpression()
    {
      return getRhsExpression_impl(std::make_index_sequence<numStoredTimes-1>{}); 
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
      return (-1.0) * m_inv_coeff * (getRhsExpression_impl_helper<idxs>() + ...); 
    }

    template<std::size_t ithNode>
    auto getRhsExpression_impl_helper()
    {
      if constexpr(std::tuple_size<MatrixTup>::value == 0){
        // just scalars need to be added 
        double scalar_sum = std::apply(
          [this](auto&&... args){
            return (args.template coeffAt<ithNode, numStoredTimes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_scalar_coeff_sum_partition
        ); 
        // need to flip the coefficient by (-1) so that its moved to rhs 
        return scalar_sum * m_stored_sols[ithNode]; 
      }
      else if constexpr(std::tuple_size<ScalarTup>::value == 0){
        // just matrix return types need to be added
        decltype(auto) matrix_sum = std::apply(
          [this](auto&&... args){
            return (args.template coeffAt<ithNode, numStoredTimes>(m_weights_calc.getArray()) + ...); 
          },
          m_mat_coeff_sum_partition
        ); 
        return matrix_sum * m_stored_sols[ithNode]; 
      }
      else{
        // sum of both scalar AND matrix return types! 
        double scalar_sum = std::apply(
          [this](auto&&... args){
            return (args.template coeffAt<ithNode, numStoredTimes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_scalar_coeff_sum_partition
        );
        decltype(auto) matrix_sum = std::apply(
          [this](auto&&... args){
            return (args.template coeffAt<ithNode, numStoredTimes>(m_weights_calc.getArray()) + ...); 
          }, 
          m_mat_coeff_sum_partition
        ); 
        return scalar_sum * m_stored_sols[ithNode] + matrix_sum * m_stored_sols[ithNode]; 
      }
    }

    // gets 1 / c where c is coeff of U(n+1) in fdm equation 
    void buildInvCoeff()
    {
      std::cout << "scalar tup size: " << std::tuple_size<ScalarTup>::value << std::endl; 
      std::cout << "matrix tup size: " << std::tuple_size<MatrixTup>::value << std::endl; 
      std::cout << "invcoeff scalar? " << std::is_same<fdm::Scalar, InvCoeff>::value << std::endl; 
      if constexpr(std::tuple_size<MatrixTup>::value == 0){
        // all coeffAt's evaluate to scalar -> return 1 / sum(coeffs) 
        fdm::Scalar s = std::apply(
            [this](auto&&... coeffs){
              return (coeffs.template coeffAt<numStoredTimes-1, numStoredTimes>(m_weights_calc.getArray()) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        m_inv_coeff = 1.0 / s;  
      }
      else if constexpr(std::tuple_size<ScalarTup>::value == 0){
        // all coeffAt's evaluate to matrix -> return (sum(coeffs)).cwiseInverse 
        auto expr = std::apply(
              [this](auto&&... coeffs){
                return (coeffs.template coeffAt<numStoredTimes-1, numStoredTimes>(m_weights_calc.getArray()) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 
        m_inv_coeff = (expr.diagonal().cwiseInverse()).asDiagonal();
      }
      else{
        // otherwise return product of 1/(sum(scalar) + (sum(Mats))
        std::cout << "correct branch hit" << std::endl; 
        double s = std::apply(
            [this](auto&&... coeffs){
              return (coeffs.template coeffAt<numStoredTimes-1, numStoredTimes>(m_weights_calc.getArray()) + ...); 
            }, 
            m_scalar_coeff_sum_partition
        );
        auto expr = std::apply(
              [this](auto&&... coeffs){
                return (coeffs.template coeffAt<numStoredTimes-1, numStoredTimes>(m_weights_calc.getArray()) + ...); 
              }, 
              m_mat_coeff_sum_partition 
        ); 
        // combine repeated s, A along diagonal 
        m_inv_coeff.diagonal() = (expr.diagonal().array() + s).inverse(); 
      }
    }

    // Sets mesh onto a Time Derivative. traverse lhs/rhs of multiply expressions
    template<typename TimeDerivU>
    void setMesh_singleton(const std::shared_ptr<const fdm::Mesh>& m,TimeDerivU& tderiv)
    {
      if constexpr(texprs::traits::is_coeffmult_crtp<TimeDerivU>::value){
        setMesh_singleton(m, tderiv.getLhs()); 
        setMesh_singleton(m, tderiv.getRhs()); 
      }
      if constexpr(!std::is_arithmetic<TimeDerivU>::value && !texprs::traits::is_timederiv_crtp<TimeDerivU>::value){
        tderiv.setMesh(m); 
      }
    }

    // Sets mesh onto a Time Derivative. traverse lhs/rhs of multiply expressions
    template<typename TimeDerivU>
    void setTime_singleton(double t,TimeDerivU& tderiv)
    {
      if constexpr(texprs::traits::is_coeffmult_crtp<TimeDerivU>::value){
        setTime_singleton(t, tderiv.getLhs()); 
        setTime_singleton(t, tderiv.getRhs()); 
      }
      if constexpr(!std::is_arithmetic<TimeDerivU>::value && !texprs::traits::is_timederiv_crtp<TimeDerivU>::value){
          tderiv.setTime(t); 
      }
    }
}; 

// entry point for giving specific # of nodes, CTAD time derivate
template<std::size_t numStoredTimes, typename TimeDeriv>
auto make_Executor(TimeDeriv& tderiv)
{
  return Executor<TimeDeriv,numStoredTimes>(tderiv); 
} 

// entry point  without nodes being specified 
template<typename TimeDeriv>
auto make_Executor(TimeDeriv& tderiv)
{
  return Executor<TimeDeriv, TimeDeriv::maxOrder+1>(tderiv); 
}

  } // end namespace texprs 
} // end namespace fdm 

#endif // LhsExecutor.hpp