// ImplicitSolver.hpp
//
// Mediator class .. 
// in order to solver PDEs with finite difference methods 
// namely implicit steps from t(n-1) to t(n) 
//
// JAF 3/4/2026 

#ifndef FDM_SOLVERS_IMPLICITSOLVER_H
#define FDM_SOLVERS_IMPLICITSOLVER_H

#include<Eigen/IterativeLinearSolvers> // BiCGSTAB sparse iterative solver  
#include "../TExprs/TExprTraits.hpp" // check LHS is time derivatives 
#include "../TExprs/Executor.hpp" // marches through time 
#include "../OutsideSteps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../OutsideSteps/OStepBase.hpp" // StepType scoped enumeration 
#include "../Utilities/RowMajorIdentityExpr.hpp"
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fdm{
  namespace solvers{ 

// template<typename LhsExpression, typename RhsExpression, typename OutsideStepsTuple>
template<
  typename LhsExpression, 
  typename RhsExpression, 
  typename OutsideStepsTuple, 
  class SparseIterativeSolver=Eigen::BiCGSTAB<fdm::CSRMatrix>
>
class ImplicitSolver
{
  public:
    // Type Defs -------------------------------------- 
    using TExpr = LhsExpression; 
    using Linop = RhsExpression; 
    
  private:
    // Member Data -------------------------------------
    LhsExpression& m_lhs; // expression of time derivatives 
    RhsExpression& m_rhs; // expression of spatial derivatives 
    using TupleCleaned = typename std::remove_reference<OutsideStepsTuple>::type; 
    TupleCleaned m_osteps; // std::tuple<> of outside steps 
    std::unique_ptr<SparseIterativeSolver> m_iterative_solver;

  public:
    // Constructors + Destructor ===========================

    ImplicitSolver()=delete; 

    ImplicitSolver(
      LhsExpression& l_init, 
      RhsExpression& r_init, 
      OutsideStepsTuple ostep_init, 
      std::unique_ptr<SparseIterativeSolver> s_init = std::make_unique<SparseIterativeSolver>()
    )
      : m_lhs(l_init), m_rhs(r_init), m_osteps(ostep_init), m_iterative_solver(std::move(s_init))
    {
      static_assert(std::is_same_v<typename SparseIterativeSolver::MatrixType, fdm::CSRMatrix>, "must use iterative solver on fmd::Matrix"); 
    }

    // not copyable!
    ImplicitSolver(const ImplicitSolver& other)=delete; 

    // moveable 
    ImplicitSolver(ImplicitSolver&& other)=default; 

    // destructor 
    ~ImplicitSolver()=default; 

    // Member Functions ===========================================================
    auto& getIterativeSolver(){ return *m_iterative_solver; }
    const auto& getIterativeSolver() const { return *m_iterative_solver; }

    template<typename M, typename C, typename Pred = LastSaver>
    auto calculate(SolverArgs<M,C> args, Pred save_policy = {}) 
    {
      // setup time context 
      auto it = std::next(args.times->cbegin(), args.initialConditions.size()-1); 
      auto time_ctx = fdm::osteps::make_time(*it, 0.0, std::move(args.times));
      ++it; 

      // setup executor 
      auto executor = fdm::texprs::make_Executor(m_lhs); 
      executor.pushSolutionRange(args.initialConditions.begin(), args.initialConditions.end()); 
      executor.pushTimeRange(time_ctx.container->cbegin(), it); 

      // set up operators. 
      fdm::utils::RowMajorIdentity identity(0,0); // will resize later  
      m_rhs.setMesh(args.mesh); 
      executor.setMesh(args.mesh); 

      // set up context 
      auto ctx = fdm::osteps::make_context(std::move(args.mesh), &executor, &m_rhs, this); 

      // store allocated memory between steps in solver hot loop  
      fdm::CSRMatrix stencil; 
      fdm::Vector rhs_vector;  
      fdm::Vector solution_u;  

      // hot loop through times
      auto end = time_ctx.container->cend();
      for(; it!= end; ++it)
      { 
        time_ctx.next = *it; 

        if constexpr(fdm::linops::internal::traits<RhsExpression>::is_timedep){
          // set the operator to the right side of the step [t(n), t(n+1)] for implicit steps 
          m_rhs.setTime(time_ctx.next); 
          // should build an autonomous solver for these linops, since it evaluate the 
          // expression at every step. but still save a little time
          // we also can't check that outside steps won't change the mesh we operate on.  
        }

        // again using left side of [t(n), t(n+1)] time step 
        executor.pushTime(time_ctx.next); 
        executor.calculate(time_ctx.next); 

        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ 
            ((lam_args.template BeforeLinAlgebra<fdm::osteps::StepType::Implicit>(time_ctx, ctx)), ...); 
          }, 
          m_osteps
        ); 
        
        // store the matrix into stencil 
        std::size_t s = m_rhs.rows(); 
        identity.resize(s,s); 
        stencil = identity - (executor.getInvCoeff() * m_rhs.toEigen());
        
        // outside steps matrix before step
        std::apply(
          [&](auto&&... lam_args){ ((lam_args.template MatBeforeStep<fdm::osteps::StepType::Implicit>(stencil, time_ctx, ctx)), ...); }, 
          m_osteps
        ); 

        // store the expression into a vector 
        rhs_vector = executor.getRhsExpression(); 

        // outside steps vector before step 
        std::apply(
          [&](auto&&... lam_args){ ((lam_args.template VecBeforeStep<fdm::osteps::StepType::Implicit>(rhs_vector, time_ctx, ctx)), ...); }, 
          m_osteps
        ); 

        // Implicit Step (I - D(t+1))*U(n+1) = rhs 
        m_iterative_solver->compute(stencil); 
        solution_u = m_iterative_solver->solveWithGuess(rhs_vector,rhs_vector);
        
        // outside steps vector after step(next_sol) 
        std::apply(
          [&](auto&&... lam_args){ ((lam_args.template VecAfterStep<fdm::osteps::StepType::Implicit>(solution_u, time_ctx, ctx)), ...); }, 
          m_osteps
        );

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
} // end namespace fdm 

#endif // ExplicitSolver.hpp 