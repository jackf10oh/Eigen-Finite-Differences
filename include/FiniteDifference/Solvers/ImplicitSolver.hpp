// ImplicitSolver.hpp
//
// Mediator class .. 
// in order to solver PDEs with finite difference methods 
// namely implicit steps from t(n-1) to t(n) 
//
// JAF 3/4/2026 

#ifndef IMPLICITSOLVER_H
#define IMPLICITSOLVER_H 

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
    std::size_t m_max_iters; // max number of iterations between time steps for iterative linear solver. 

  public:
    // Constructors + Destructor ===========================

    ImplicitSolver()=delete; 

    ImplicitSolver(LhsExpression& l_init, RhsExpression& r_init, OutsideStepsTuple ostep_init)
      : m_lhs(l_init), m_rhs(r_init), m_osteps(ostep_init), m_max_iters(20)
    {}

    // not copyable!
    ImplicitSolver(const ImplicitSolver& other)=delete; 

    // moveable 
    ImplicitSolver(ImplicitSolver&& other)=default; 

    // destructor 
    ~ImplicitSolver()=default; 

    // Member Functions ===========================================================
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
      Eigen::VectorXd rhs_vector;  
      Eigen::VectorXd solution_u;  

      // hot loop through times
      auto end = time_ctx.container->cend();
      for(; it!= end; ++it)
      { 
        time_ctx.next = *it; 

        if constexpr(m_rhs.isTimeDep){
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
        decltype(auto) xpr = m_rhs.asMatrix(); // could be expensive if rhs just gives a matrix
        std::size_t s = xpr.rows(); 
        identity.resize(s,s); 
        stencil = identity.asMatrix() - executor.getInvCoeff() * xpr;
        
        // outside steps matrix before step
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<fdm::osteps::StepType::Implicit>(stencil, time_ctx, ctx)), ...); }, 
          m_osteps
        ); 

        // store the expression into a vector 
        rhs_vector = executor.getRhsExpression(); 

        // outside steps vector before step 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template VecBeforeStep<fdm::osteps::StepType::Implicit>(rhs_vector, time_ctx, ctx)), ...); }, 
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

  } // end namespace solvers
} // end namespace fdm 

#endif // ExplicitSolver.hpp 