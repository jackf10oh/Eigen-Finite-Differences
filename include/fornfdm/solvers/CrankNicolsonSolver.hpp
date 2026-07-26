// CrankNicolsonSolver.hpp
//
// implements Crank Nicolson scheme
// Dt(U) = 1/2 [ Dx( u(n), t(n)) + Dx( u(n+1), t(n+1))] 
//
// JAF 4/12/2026 

#ifndef FORNFDM_SOLVERS_CRANKNICOLSONSOLVER_H
#define FORNFDM_SOLVERS_CRANKNICOLSONSOLVER_H

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
class CrankNicolsonSolver : public SolverBase<CrankNicolsonSolver<LhsType, RhsType, OStepTup, SparseIterativeSolver>, LhsType, RhsType, OStepTup>
{
  public:
    // Type Defs -------------------------------------- 
    using Base = SolverBase<CrankNicolsonSolver<LhsType, RhsType, OStepTup,SparseIterativeSolver>, LhsType, RhsType, OStepTup>;

  private:
    // Member Data -------------------------------------
    std::unique_ptr<SparseIterativeSolver> m_iterative_solver;

  public:
    // Constructors + Destructor ==============================================================

    CrankNicolsonSolver()=delete; 

    CrankNicolsonSolver(
      LhsType& l_init, 
      RhsType& r_init, 
      OStepTup ostep_init, 
      std::unique_ptr<SparseIterativeSolver> s_init = std::make_unique<SparseIterativeSolver>()
    )
      : Base(l_init, r_init, std::move(ostep_init)), m_iterative_solver(std::move(s_init))
    {
      static_assert(std::is_same_v<typename SparseIterativeSolver::MatrixType, fornfdm::CSRMatrix>, "must use iterative solver on fmd::Matrix"); 
    }

    // not copyable! 
    CrankNicolsonSolver(const CrankNicolsonSolver& other)=delete; 

    // moveable 
    CrankNicolsonSolver(CrankNicolsonSolver&& other)=default; 

    // destructor 
    ~CrankNicolsonSolver()=default; 

    // Member Functions ==============================================================
    auto& getIterativeSolver(){ return *m_iterative_solver; }
    const auto& getIterativeSolver() const { return *m_iterative_solver; }

    template<typename M, typename C, typename Pred = LastSaver>
    auto calculate(SolverArgs<M,C> args, Pred save_policy = {})  
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
      fornfdm::CSRMatrix linop_untouched; 
      fornfdm::Vector rhs_vector;  
      fornfdm::Vector solution_u;  

      // need to get first linop at time.now before it caches the next time steps... 
      linop_untouched = 0.5 * (this->m_rhs.evalTime(time_ctx.now)); 

      // hot loop through times
      auto end = time_ctx.container->cend();
      for(; it!= end; ++it)
      { 
        time_ctx.next = *it; 

        // again using left side of [t(n), t(n+1)] time step 
        executor.pushTime(time_ctx.next); 
        executor.calculate(0.5 * (time_ctx.now + time_ctx.next)); 

        // outside steps before any type of linear algebra is performed... 
        this->template tupleBeforeLinAlgebra<fornfdm::osteps::StepType::Implicit>(time_ctx,ctx);

        // store the expression into a vector 
        rhs_vector = (executor.getInvCoeff() * (linop_untouched * executor.getCurrentSolution())) + executor.getRhsExpression(); 

        // outside steps vector before step 
        this->template tupleApplyBeforeVec<fornfdm::osteps::StepType::Implicit>(rhs_vector, time_ctx, ctx);  
        
        // store the matrix into stencil 
        linop_untouched = 0.5 * (this->m_rhs.evalTime(time_ctx.next)); 
        std::size_t s = linop_untouched.rows(); 
        stencil = fornfdm::utils::Identity(s,s) - executor.getInvCoeff() * linop_untouched; 
        
        // outside steps matrix before step
        this->template tupleApplyBeforeMat<fornfdm::osteps::StepType::Implicit>(stencil, time_ctx, ctx); 

        // Implicit Step (I - D(t+1))*U(n+1) = rhs 
        m_iterative_solver->compute(stencil); 
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

template<typename L, typename R, typename O, class S=Eigen::BiCGSTAB<fornfdm::CSRMatrix>>
using CNSolver = CrankNicolsonSolver<L,R,O,S>; 

  } // end namespace solvers
} // end namespace fornfdm 

#endif // CrankNicolsonSolver.hpp 