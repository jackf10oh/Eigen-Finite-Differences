// ImplicitSolver.hpp
//
// Mediator class .. 
// in order to solver PDEs with finite difference methods 
// namely implicit steps from t(n-1) to t(n) 
//
// JAF 3/4/2026 

#ifndef IMPLICITSOLVER_H
#define IMPLICITSOLVER_H 

#include "../LinOps/LinOpTraits.hpp" // check RHS is 1D or XD LinOp + fdm::Matrix 
#include "../TExprs/TExprTraits.hpp" // check LHS is time derivatives 
#include "../TExprs/Executor.hpp" // marches through time 
#include "../OutsideSteps/StepContexts.hpp"  // feed to outside steps tuple 
#include "../OutsideSteps/OStepBase.hpp" // StepType scoped enumeration 
#include "SolverArgs.hpp"
#include "WritePolicies.hpp"

namespace fdm{
  namespace solvers{ 

template<typename LHS_EXPR, typename RHS_EXPR, typename OSTEP_TUP>
class ImplicitSolver
{
  private:
    // Member Data -------------------------------------
    LHS_EXPR& m_lhs; // expression of time derivatives 
    RHS_EXPR& m_rhs; // expression of spatial derivatives 
    typename std::remove_reference<OSTEP_TUP>::type m_ostep_tup; 
    std::size_t m_max_iters; // max number of iterations between time steps for iterative linear solver. 
  public:
    // Constructors + Destructor ===========================
    ImplicitSolver()=delete; 
    ImplicitSolver(LHS_EXPR& l_init, RHS_EXPR& r_init, OSTEP_TUP ostep_init)
      : m_lhs(l_init), m_rhs(r_init), m_ostep_tup(ostep_init), m_max_iters(20)
    {}
    ImplicitSolver(const ImplicitSolver& other)=delete; 
    ~ImplicitSolver()=default; 

    // Member Functions 
    void SetMaxIterations(std::size_t i){ m_max_iters = i; } 
    void MaxIterations(std::size_t i){ m_max_iters = i; } 

    template<typename M, typename WRITE_POLICY_T = LastSaver, template<typename MAT_T> class EIGENSOLVER_T=Eigen::BiCGSTAB>
    auto Calculate(SolverArgs<M> args, WRITE_POLICY_T save_policy = {}) const 
    {
      texprs::TExprExecutor exec(m_lhs); 
      EIGENSOLVER_T<texprs::MatrixStorage_t> iterative_solver; // Eigen sparse iterative solver
      iterative_solver.setMaxIterations(m_max_iters); 
      texprs::MatrixStorage_t Mat; 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = std::next(args.time_mesh_ptr->cbegin(), args.ICs.size() - 1); 
      auto end = std::prev(args.time_mesh_ptr->cend());

      exec.ConsumeSolutionList(args.ICs.cbegin(), args.ICs.cend()); 
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), std::next(it)); 

      double t = *it; 
      while(it!= end)
      {
        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ ((lam_args.template BeforeLinAlgebra<decltype(it),decltype(exec),decltype(m_rhs),OSteps::StepType::Implicit>(it, args.domain_mesh_ptr, exec, m_rhs)), ...); }, 
          m_ostep_tup
        ); 

        // working on right end of [t(n-1), t(n)]
        ++it;
        t = *it; 

        // moved into an ostep that sets time + mesh of lhs executor / rhs expression 
        // m_rhs.setTime(t);
        // exec.setTime(t); 

        exec.BuildNextTime(t); 
        
        // Mat = I - inv_coeff * FDStencil ; 
        Mat = m_rhs.getMat(); 
        // scale Mat according to 1 / dt ... 
        if constexpr(std::is_same<decltype(exec.inv_coeff()), const double&>::value){
          // INV_COEFF_T is a scalar
          Mat  *= -exec.inv_coeff(); 
        }
        else{
          // INV_COEFF_T is a Matrix 
          texprs::MatrixStorage_t temp = -exec.inv_coeff() * Mat; 
          Mat = std::move(temp); 
        }
        for(std::size_t i=0; i<Mat.rows(); ++i) Mat.coeffRef(i,i) += 1.0; 

        // outside steps matrix before step(Mat) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<OSteps::StepType::Implicit>(t, args.domain_mesh_ptr, Mat)), ...); }, 
          m_ostep_tup
        ); 

        Eigen::VectorXd rhs = std::move(exec.RhsVector()); 
        // outside steps solution before step (rhs) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolBeforeStep<OSteps::StepType::Implicit>(t, args.domain_mesh_ptr, rhs)), ...); }, 
          m_ostep_tup
        ); 

        // Explicit Step ( U(n+1) = D*U(n) + U(n) )
        // Implicit Step (I - D(t+1))*U(n+1) = rhs 
        iterative_solver.compute(Mat); 
        Eigen::VectorXd next_sol = iterative_solver.solveWithGuess(rhs,rhs); 
        
        // ourside steps solution after step(next_sol) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolAfterStep<OSteps::StepType::Implicit>(t, args.domain_mesh_ptr, next_sol)), ...); }, 
          m_ostep_tup
        ); 

        // give expiring solution to WRITE_POLICY_T
        save_policy.saveSolution(std::move(exec.ExpiringSol())); 

        // push Solution, time to executor 
        exec.ConsumeSolution(next_sol);
        exec.ConsumeTime(t); 
      }    

      // Write remaining solutions to save_policy 
      for(auto i=0; i<exec.StoredSols().size()-1; i++) save_policy.saveSolution( std::move(exec.StoredSols()[i])); 

      // write policy also determines return type / how to handle last solution
      return save_policy.saveLastSolution(std::move( exec.MostRecentSol() )); 
    } // end .CalculateImp(args, save_policy) 

}; 

  } // end namespace solvers
} // end namespace fdm 

#endif // ExplicitSolver.hpp 