// FastExpSolver.hpp
//
// class that ducks some of the main contracts in Solvers interface. 
// notably there is no beforeLinAlgebra, applyBeforeMat, or applyBeforeVec
// the full CSRMatrix is never calculated, so RHS operator is possibly time dependent.
//
// JAF 8/4/2026

#ifndef FORNFDM_SOLVERS_FastExpSolver_H
#define FORNFDM_SOLVERS_FastExpSolver_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include "../diffops/traits.hpp" // check RHS is_timedep
#include "../texprs/traits.hpp" // check LHS is time derivatives 
#include "../texprs/Executor.hpp" // marches through time 
#include "../outside_steps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../outside_steps/OStepBase.hpp" // StepType scoped enumeration 
#include "SolverBase.hpp"
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fornfdm{
  namespace solvers{ 

template<typename LhsType, typename RhsType, typename OStepTup>
class FastExpSolver : public SolverBase<FastExpSolver<LhsType, RhsType, OStepTup>, LhsType, RhsType, OStepTup> 
{
  private:
    // Type Defs -------
    using Base = SolverBase<FastExpSolver<LhsType, RhsType, OStepTup>, LhsType, RhsType, OStepTup>; 
  
    public:
    // Constructors + Destructor ==============================================

    FastExpSolver()=delete; 

    FastExpSolver(LhsType& l_init, RhsType& r_init, OStepTup ostep_init)
      : Base(l_init, r_init, std::move(ostep_init))
    {}

    // not copyable! 
    FastExpSolver(const FastExpSolver& other)=delete; 

    // moveable
    FastExpSolver(FastExpSolver&& other)=default;  

    // destructor  
    ~FastExpSolver()=default; 

    // Member Functions ======================================================
    template<class M, class C, class Pred = LastSaver>
    auto calculate(SolverArgs<M, C> args, Pred save_policy = {}) 
    {
      // Fast Solvers only accept a SharedConstMesh and const TimeArgs for arguments!
      static_assert(std::is_same<M, const fornfdm::Mesh>::value, "FastImpSolver only takes const Mesh for domain.");
      static_assert(std::is_same<C, const solvers::TimeArg>::value, "FastImpSolver only takes TimeArg for times.");

      fornfdm::Real t_start = args.times->getStart();
      fornfdm::Real t_stop = args.times->getStop();
      fornfdm::Real t_stepsize = args.times->getStepSize();
      std::size_t end = args.times->getNumSteps();


      // setup time context 
      std::size_t step_n = args.initialConditions.size()-1;
      auto time_ctx = fornfdm::osteps::make_time(t_start + step_n*t_stepsize, 0.0, std::move(args.times));
      ++step_n; 

      // setup executor 
      auto executor = fornfdm::texprs::make_Executor(this->m_lhs); 

      auto n = step_n - executor.numStoredTimes + 1;

      for(auto i=1; i<executor.numStoredTimes; ++i){
        executor.getStoredTimes()[i] = t_start + t_stepsize * n;
        ++n;
      };

      auto sol_end = executor.pushSolutionRange(args.initialConditions.begin(), args.initialConditions.end()); 
      for(auto sol_it=args.initialConditions.begin(); sol_it != sol_end; ++sol_it)
      {
        save_policy.saveSolution(std::move(*sol_it)); 
      }

      // set up operators. 
      this->m_rhs.setMesh(args.mesh); 
      executor.setMesh(args.mesh); 
      executor.calculate(t_start + step_n * t_stepsize);

      // set up context 
      auto ctx = fornfdm::osteps::make_context(std::move(args.mesh), &executor, &(this->m_rhs), this); 

      fornfdm::Vector solution_u;  
      // hot loop through times
      for(; step_n != end; ++step_n)
      { 
        time_ctx.next = t_start + step_n * t_stepsize; 
        executor.pushTime(time_ctx.next);
        executor.calculate(time_ctx.next);
        this->m_rhs.setTime(time_ctx.now);

        // no beforeLinAlgebra!
        // no applyBeforeMat!
        // no applyBeforeVec! 

        // Explicit Step 
        solution_u = (executor.getInvCoeff() * this->m_rhs.toEigen()) * executor.getCurrentSolution() + executor.getRhsExpression(); 
        
        // outside steps solution after step(next_sol) 
        this->template tupleApplyAfterVec<fornfdm::osteps::StepType::Explicit>(solution_u, time_ctx, ctx); 

        // save oldest solution before it goes 
        save_policy.saveSolution(executor.getExpiringSolution()); 

        // push solution into executor
        executor.rotateStoredSolutions(1); 
        executor.getStoredSolutions().back().swap(solution_u); 

        // advance time_ctx 1 step in time 
        time_ctx.now = time_ctx.next; 
      }        
      // Give remaining solutions to save_policy 
      for(auto i=0; i < executor.numStoredSols-1; ++i) save_policy.saveSolution( std::move(executor.getStoredSolutions()[i])); 

      // save_policy also determines return type / how to handle last solution
      return save_policy.saveLastSolution(std::move( executor.getCurrentSolution() )); 
    } 
}; 

  } // end namespace solvers
} // end namespace fornfdm 

#endif // FastExpSolver.hpp 