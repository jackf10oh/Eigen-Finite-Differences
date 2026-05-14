// SolverBase.hpp
//
// Base class between all solvers
// to give a uniform interface 
// and some protected utilities to derived 
//
// JAF 5/12/2026 

#ifndef FORNFDM_SOLVERS_SOLVERBASE_H
#define FORNFDM_SOLVERS_SOLVERBASE_H 

#include"../Types.hpp" // Matrix, Vector
#include "../OutsideSteps/OStepBase.hpp" // StepType 
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fornfdm{
  namespace solvers{

template<typename Derived, typename LhsType, typename RhsType, typename OStepTup>
class SolverBase
{
  public:
    // Member Data -------------------------------------

    /* minimum is always fixed == to order of Lhs time deriv expression. 
    higher order schemes. i.e. Runge-Kutta, Adams-Bashforth
    will have to accept this minimum and use a different scheme to fill 
    to their respective minimum */
    static constexpr std::size_t numIcsMin = LhsType::maxOrder;  
  
  protected:
    LhsType& m_lhs; // expression of time derivatives 
    RhsType& m_rhs; // expression of spatial derivatives 
    using TupleCleaned = typename std::remove_reference<OStepTup>::type; 
    TupleCleaned m_osteps; // std::tuple<> of outside steps 

  public:
    // Constructors + Destructor ==============================================

    SolverBase()=delete; 

    SolverBase(LhsType& l_init, RhsType& r_init, OStepTup ostep_init)
      : m_lhs(l_init), m_rhs(r_init), m_osteps(std::move(ostep_init))
    {}

    // not copyable! 
    SolverBase(const SolverBase& other)=delete; 

    // moveable
    SolverBase(SolverBase&& other)=default;  

    // destructor  
    ~SolverBase()=default; 

    // Member Functions ======================================================
    auto& derived(){ return *static_cast<Derived*>(this); }
    const auto& derived() const { return *static_cast<const Derived*>(this); }
    template<class M, class C, class Pred = LastSaver>
    auto calculate(SolverArgs<M, C> args, Pred save_policy = {})
    { 
      return derived().calculate(std::move(args),save_policy); 
    }

  protected:
    template<fornfdm::osteps::StepType step, class TCtx, class Ctx>
    void tupleBeforeLinAlgebra(const TCtx& time_ctx, Ctx& ctx)
    {
        std::apply(
          [&](auto&&... lam_args){ 
            ((lam_args.template BeforeLinAlgebra<step>(time_ctx, ctx)), ...); 
          }, 
          m_osteps
        ); 
    }

    template<fornfdm::osteps::StepType step, class TCtx, class Ctx>
    void tupleMatBeforeStep(fornfdm::CSRMatrix& mat, const TCtx& time_ctx, const Ctx& ctx)
    {
      std::apply(
        [&](auto&&... lam_args)
        { 
          ((lam_args.template MatBeforeStep<step>(mat, time_ctx, ctx)), ...); 
        }, 
        m_osteps
      ); 
    }

    template<fornfdm::osteps::StepType step, class TCtx, class Ctx>
    void tupleVecBeforeStep(fornfdm::Vector& rhs_vector, const TCtx& time_ctx, const Ctx& ctx)
    {
      // outside steps vector before step 
      std::apply(
        [&](auto&&... lam_args){ ((lam_args.template VecBeforeStep<fornfdm::osteps::StepType::Explicit>(rhs_vector, time_ctx, ctx)), ...); }, 
        m_osteps
      ); 
    }

    template<fornfdm::osteps::StepType step, class TCtx, class Ctx>
    void tupleVecAfterStep(fornfdm::Vector solution_u, const TCtx& time_ctx, const Ctx& ctx)
    {
      // outside steps solution after step(next_sol) 
      std::apply(
        [&](auto&... lam_args)
        { 
          ((lam_args.template VecAfterStep<step>(solution_u, time_ctx, ctx)), ...); 
        }, 
        m_osteps
      );
    }
}; 

  } // end namespace solvers
} // end namespace fornfdm

#endif // SolverBase.hpp 