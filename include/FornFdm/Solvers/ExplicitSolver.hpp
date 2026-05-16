// ExplicitSolver.hpp
//
// Mediator class that manages 3 classes:
// 1. An expression of time derivatives (TExpr) 
// 2. An expression of spatial derivatives (LinOps)
// 3. A tuple of Outside steps (OSteps) 
// in order to solver PDEs with finite difference methods 
// namely explicit steps from t(n-1) to t(n) 
//
// JAF 3/4/2026 

#ifndef FORNFDM_SOLVERS_EXPLICITSOLVER_H
#define FORNFDM_SOLVERS_EXPLICITSOLVER_H

#include "../TExprs/Traits.hpp" // check LHS is time derivatives 
#include "../TExprs/Executor.hpp" // marches through time 
#include "../OutsideSteps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../OutsideSteps/OStepBase.hpp" // StepType scoped enumeration 
#include "SolverBase.hpp"
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fornfdm{
  namespace solvers{ 

template<typename LhsType, typename RhsType, typename OStepTup>
class ExplicitSolver : public SolverBase<ExplicitSolver<LhsType, RhsType, OStepTup>, LhsType, RhsType, OStepTup> 
{
  private:
    // Type Defs -------
    using Base = SolverBase<ExplicitSolver<LhsType, RhsType, OStepTup>, LhsType, RhsType, OStepTup>; 
  
    public:
    // Constructors + Destructor ==============================================

    ExplicitSolver()=delete; 

    ExplicitSolver(LhsType& l_init, RhsType& r_init, OStepTup ostep_init)
      : Base(l_init, r_init, std::move(ostep_init))
    {}

    // not copyable! 
    ExplicitSolver(const ExplicitSolver& other)=delete; 

    // moveable
    ExplicitSolver(ExplicitSolver&& other)=default;  

    // destructor  
    ~ExplicitSolver()=default; 

    // Member Functions ======================================================
    template<class M, class C, class Pred = LastSaver>
    auto calculate(SolverArgs<M, C> args, Pred save_policy = {}) 
    {
      // setup time context 
      auto it = std::next(args.times->cbegin(), args.initialConditions.size()-1); 
      auto time_ctx = fornfdm::osteps::make_time(*it, 0.0, std::move(args.times));
      ++it; 

      // setup executor 
      auto executor = fornfdm::texprs::make_Executor(this->m_lhs); 
      executor.pushTimeRange(time_ctx.container->cbegin(), it); 
      auto sol_end = executor.pushSolutionRange(args.initialConditions.begin(), args.initialConditions.end()); 
      for(auto sol_it=args.initialConditions.begin(); sol_it != sol_end; ++sol_it)
      {
        // Eventually want save policies to take by (solution, time) 
        save_policy.saveSolution(std::move(*sol_it)); 
      }

      // set up operators. 
      this->m_rhs.setMesh(args.mesh); 
      executor.setMesh(args.mesh); 

      // set up context 
      auto ctx = fornfdm::osteps::make_context(std::move(args.mesh), &executor, &(this->m_rhs), this); 

      // store allocated memory between steps in solver hot loop  
      fornfdm::CSRMatrix stencil; 
      fornfdm::Vector rhs_vector;  
      fornfdm::Vector solution_u;  

      // hot loop through times
      auto end = time_ctx.container->cend();
      for(; it!= end; ++it)
      { 
        time_ctx.next = *it; 

        if constexpr(fornfdm::linops::internal::traits<RhsType>::is_timedep){
          // set the operator to the left side of the step [t(n), t(n+1)] for explicit steps 
          this->m_rhs.setTime(time_ctx.now); 
          // should build an autonomous solver for these linops, since it evaluate the 
          // expression at every step. but still save a little time
          // we also can't check that outside steps won't change the mesh we operate on.  
        }

        // again using left side of [t(n), t(n+1)] time step 
        executor.pushTime(time_ctx.next); 
        executor.calculate(time_ctx.now); 

        // outside steps before any type of linear algebra is performed... 
        this->template tupleBeforeLinAlgebra<fornfdm::osteps::StepType::Explicit>(time_ctx,ctx); 
        
        // store the matrix into stencil 
        stencil = executor.getInvCoeff() * this->m_rhs.toEigen();
        
        // outside steps matrix before step
        this->template tupleMatBeforeStep<fornfdm::osteps::StepType::Explicit>(stencil, time_ctx, ctx); 

        // store the expression into a vector 
        rhs_vector = executor.getRhsExpression(); 

        // outside steps vector before step 
        this->template tupleVecBeforeStep<fornfdm::osteps::StepType::Explicit>(rhs_vector, time_ctx, ctx); 

        // Explicit Step 
        solution_u = stencil * executor.getCurrentSolution() + rhs_vector; 
        
        // outside steps solution after step(next_sol) 
        this->template tupleVecAfterStep<fornfdm::osteps::StepType::Explicit>(solution_u, time_ctx, ctx); 

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

#endif // ExplicitSolver.hpp 