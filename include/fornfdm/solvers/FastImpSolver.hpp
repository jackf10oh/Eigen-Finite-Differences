// FastExpSolver.hpp
//
// Class that ducks some of the main contracts in Solvers interface. 
// Notably, there is no beforeLinAlgebra, or applyBeforeMat 
// The full CSRMatrix is never calculated once upfront. 
// So RHS operator MUST be autonomous.
//
// JAF 8/4/2026


#ifndef FORNFDM_SOLVERS_FASTIMPSOLVER_H
#define FORNFDM_SOLVERS_FASTIMPSOLVER_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include<Eigen/IterativeLinearSolvers> // BiCGSTAB sparse iterative solver  
#include "../diffops/traits.hpp" // check RHS is_timedep
#include "../texprs/traits.hpp" // check LHS is time derivatives 
#include "../texprs/Executor.hpp" // marches through time 
#include "../outside_steps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../outside_steps/OStepBase.hpp" // StepType scoped enumeration 
#include "../utilities/Identity.hpp"
#include "SolverBase.hpp"
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fornfdm{
  namespace solvers{ 

template<
  typename LhsType, 
  typename RhsType, 
  typename OStepTup, 
  class SparseIterativeSolver=Eigen::BiCGSTAB<fornfdm::CSRMatrix>
>
class FastImpSolver : public SolverBase<FastImpSolver<LhsType, RhsType, OStepTup,SparseIterativeSolver>, LhsType, RhsType, OStepTup> 
{
  private:
    // Type Defs -------------------------------------- 
    using Base = SolverBase<FastImpSolver<LhsType, RhsType, OStepTup,SparseIterativeSolver>, LhsType, RhsType, OStepTup>;   
    
  private:
    // Member Data -------------------------------------
    std::unique_ptr<SparseIterativeSolver> m_iterative_solver;

  public:
    // Constructors + Destructor ===========================

    FastImpSolver()=delete; 

    FastImpSolver(
      LhsType& l_init, 
      RhsType& r_init, 
      OStepTup ostep_init, 
      std::unique_ptr<SparseIterativeSolver> s_init = std::make_unique<SparseIterativeSolver>()
    )
      : Base(l_init, r_init, std::move(ostep_init)), m_iterative_solver(std::move(s_init))
    {
      static_assert(std::is_same_v<typename SparseIterativeSolver::MatrixType, fornfdm::CSRMatrix>, "must use iterative solver on fornfdm::CSRMatrix"); 
      static_assert(fornfdm::linops::internal::traits<RhsType>::is_timedep == false, "error. FastImpSolver calculates stencil only once upfront. RHS must be autonomous!");
    }

    // not copyable!
    FastImpSolver(const FastImpSolver& other)=delete; 

    // moveable 
    FastImpSolver(FastImpSolver&& other)=default; 

    // destructor 
    ~FastImpSolver()=default; 

    // Member Functions ===========================================================
    auto& getIterativeSolver(){ return *m_iterative_solver; }
    const auto& getIterativeSolver() const { return *m_iterative_solver; }

    template<class M, class C, class Pred = LastSaver>
    auto calculate(SolverArgs<M, C> args, Pred save_policy = {}) 
    {
      // Fast Solvers only accept a SharedConstMesh and const TimeArgs for arguments!
      static_assert(std::is_same<M, const fornfdm::Mesh>::value, "FastExpSolver only takes const Mesh for domain.");
      static_assert(std::is_same<C, const solvers::TimeArg>::value, "FastExpSolver only takes TimeArg for times.");

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

      for(auto i=0; i<executor.numStoredTimes; ++i){
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
      executor.calculate(t_start + t_stepsize * step_n);

      // Build the stencil. 
      std::size_t s = args.mesh->sizesProduct();
      fornfdm::CSRMatrix stencil = fornfdm::utils::Identity(s,s) - executor.getInvCoeff() * this->m_rhs.toEigen();

      // set up context 
      auto ctx = fornfdm::osteps::make_context(std::move(args.mesh), &executor, &(this->m_rhs), this); 

      // applyBeforeMat only fires once!
      this->template tupleApplyBeforeMat<fornfdm::osteps::StepType::Implicit>(stencil, time_ctx, ctx);

      // don't need to calculate preconditioner each time.
      m_iterative_solver->compute(stencil); 

      fornfdm::Vector rhs_vector;
      fornfdm::Vector solution_u;

      // hot loop through times
      for(; step_n != end; ++step_n)
      { 
        time_ctx.next = t_start + step_n * t_stepsize; 

        // no beforeLinAlgebra!
        // no applyBeforeMat!

        rhs_vector = executor.getRhsExpression();
        
        // applyBeforeVec
        this->template tupleApplyBeforeVec<fornfdm::osteps::StepType::Implicit>(rhs_vector, time_ctx, ctx);

        // Implicit Step 
        solution_u = m_iterative_solver->solveWithGuess(rhs_vector,rhs_vector); 
        
        // outside steps solution after step(next_sol) 
        this->template tupleApplyAfterVec<fornfdm::osteps::StepType::Implicit>(solution_u, time_ctx, ctx); 

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

#endif // FastImpSolver.hpp 