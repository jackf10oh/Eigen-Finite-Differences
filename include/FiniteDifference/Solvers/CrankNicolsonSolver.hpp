// CrankNicolsonSolver.hpp
//
// implements Crank Nicolson scheme
// Dt(U) = 1/2 [ Dx( u(n), t(n)) + Dx( u(n+1), t(n+1))] 
//
// JAF 4/12/2026 

#ifndef CRANKNICOLSONSOLVER_H
#define CRANKNICOLSONSOLVER_H 

#include<Eigen/IterativeLinearSolvers> // BiCGSTAB sparse iterative solver  
#include "../LinOps/LinOpTraits.hpp" // check RHS is 1D or XD LinOp + fdm::Matrix 
#include "../TExprs/TExprTraits.hpp" // check LHS is time derivatives 
#include "../TExprs/Executor.hpp" // marches through time 
#include "../OutsideSteps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../OutsideSteps/OStepBase.hpp" // StepType scoped enumeration 
#include "SolverArgs.hpp"
#include "SavePolicies.hpp"

namespace fdm{
  namespace solvers{ 

// template<typename LhsExpression, typename RhsExpression, typename OutsideStepsTuple>
template<typename LhsExpression, typename RhsExpression, typename OutsideStepsTuple>
class CrankNicolsonSolver
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
    std::size_t m_max_iters; // max number of iterations between time steps for iterative linear solver. 

  public:
    // Constructors + Destructor ==============================================================

    CrankNicolsonSolver()=delete; 

    CrankNicolsonSolver(LhsExpression& l_init, RhsExpression& r_init, OutsideStepsTuple ostep_init)
      : m_lhs(l_init), m_rhs(r_init), m_osteps(ostep_init), m_max_iters(20)
    {}

    // not copyable! 
    CrankNicolsonSolver(const CrankNicolsonSolver& other)=delete; 

    // moveable 
    CrankNicolsonSolver(CrankNicolsonSolver&& other)=default; 

    // destructor 
    ~CrankNicolsonSolver()=default; 

    // Member Functions ==============================================================
    void setMaxIterations(std::size_t i){ m_max_iters = i; } 
    auto getMaxIterations() const { return m_max_iters; } 

    template<typename M, typename C, typename Pred = LastSaver, template<typename Matrix> class SparseIterativeSolver=Eigen::BiCGSTAB>
    auto calculate(
      SolverArgs<M,C> args, 
      Pred save_policy = {}, 
      SparseIterativeSolver<fdm::Matrix> iterative_solver = {}
    ) const 
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
      fdm::linops::IOp identity; 
      m_rhs.setMesh(args.mesh); 
      executor.setMesh(args.mesh); 

      // set up context 
      auto ctx = fdm::osteps::make_context(std::move(args.mesh), &executor, &m_rhs, this); 

      // store allocated memory between steps in solver hot loop  
      fdm::Matrix stencil; 
      fdm::Matrix linop_untouched; 
      Eigen::VectorXd rhs_vector;  
      Eigen::VectorXd solution_u;  

      // need to get first linop at time.now before it caches the next time steps... 
      if constexpr(m_rhs.isTimeDep){
        m_rhs.setTime(time_ctx.now); 
      }
      linop_untouched = 0.5 * m_rhs.asMatrix(); 

      // hot loop through times
      auto end = time_ctx.container->cend();
      for(; it!= end; ++it)
      { 
        time_ctx.next = *it; 

        // again using left side of [t(n), t(n+1)] time step 
        executor.pushTime(time_ctx.next); 
        executor.calculate(0.5 * (time_ctx.now + time_ctx.next)); 

        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ 
            ((lam_args.template BeforeLinAlgebra<fdm::osteps::StepType::Implicit>(time_ctx, ctx)), ...); 
          }, 
          m_osteps
        ); 

        // store the expression into a vector 
        rhs_vector = (executor.getInvCoeff() * (linop_untouched * executor.getCurrentSolution())) + executor.getRhsExpression(); 

        // outside steps vector before step 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template VecBeforeStep<fdm::osteps::StepType::Implicit>(rhs_vector, time_ctx, ctx)), ...); }, 
          m_osteps
        ); 
        
        // store the matrix into stencil 
        if constexpr(m_rhs.isTimeDep){
          // set the operator to the right side of the step [t(n), t(n+1)] for implicit steps 
          m_rhs.setTime(time_ctx.next); 
          // should build an autonomous solver for these linops, since it evaluate the 
          // expression at every step. but still save a little time
          // we also can't check that outside steps won't change the mesh we operate on.  
        }
        linop_untouched = 0.5 * m_rhs.asMatrix(); 
        std::size_t s = linop_untouched.rows(); 
        identity.resize(s,s); 
        stencil = identity.asMatrix() - executor.getInvCoeff() * linop_untouched; 
        
        // outside steps matrix before step
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<fdm::osteps::StepType::Implicit>(stencil, time_ctx, ctx)), ...); }, 
          m_osteps
        ); 

        // Implicit Step (I - D(t+1))*U(n+1) = rhs 
        iterative_solver.setMaxIterations(m_max_iters); 
        iterative_solver.compute(stencil); 
        Eigen::VectorXd solution_u = iterative_solver.solveWithGuess(rhs_vector,rhs_vector);
        
        // outside steps vector after step(next_sol) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template VecAfterStep<fdm::osteps::StepType::Implicit>(solution_u, time_ctx, ctx)), ...); }, 
          m_osteps
        );

        // save oldest solution before it goes 
        save_policy.saveSolution(executor.getExpiringSolution()); 

        // push solution into executor
        executor.pushSolution(std::move(solution_u));

        // advance time_ctx 1 step in time 
        time_ctx.now = time_ctx.next; 
      }        
      // Give remaining solutions to save_policy 
      for(auto i=0; i < executor.numStoredSols-1; ++i) save_policy.saveSolution( std::move(executor.getStoredSolutions()[i])); 

      // save_policy also determines return type / how to handle last solution
      return save_policy.saveLastSolution(std::move( executor.getCurrentSolution() )); 
    } 
  
}; 

template<typename L, typename R, typename O>
using CNSolver = CrankNicolsonSolver<L,R,O>; 

  } // end namespace solvers
} // end namespace fdm 

#endif // CrankNicolsonSolver.hpp 