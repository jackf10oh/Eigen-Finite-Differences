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

#ifndef EXPLICITSOLVER_H
#define EXPLICITSOLVER_H 

#include<LinOps/LinOpTraits.hpp> // trait to check RHS has .GetMat() 
#include<TExprs/TimeDerivBase.hpp> // trait to check LHS is a time deriv + MatrixStorage_t 
#include<TExprs/TExprExecutor.hpp> // TExprExecutor class 
#include<OutsideSteps/OStepBase.hpp> // FDStep_Type scoped enumeration 
#include "SolverArgs.hpp"
#include "WritePolicies.hpp"

namespace Solvers{ 

template<typename LHS_EXPR, typename RHS_EXPR, typename OSTEP_TUP>
class ExplicitSolver
{
  private:
    // Member Data -------------------------------------
    LHS_EXPR& m_lhs; // expression of time derivatives 
    RHS_EXPR& m_rhs; // expression of spatial derivatives 
    typename std::remove_reference<OSTEP_TUP>::type m_ostep_tup; 

  public:
    // Constructors + Destructor ===========================
    ExplicitSolver()=delete; 
    ExplicitSolver(LHS_EXPR& l_init, RHS_EXPR& r_init, OSTEP_TUP ostep_init)
      : m_lhs(l_init), m_rhs(r_init), m_ostep_tup(ostep_init)
    {}
    ExplicitSolver(const ExplicitSolver& other)=delete; 
    ~ExplicitSolver()=default; 

    // Member Functions 
    template<typename M, typename WRITE_POLICY_T = FinalWrite>
    auto Calculate(SolverArgs<M> args, WRITE_POLICY_T write_policy = {}) const 
    {

      TExprs::TExprExecutor exec(m_lhs); 
      TExprs::MatrixStorage_t Mat; 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = std::next(args.time_mesh_ptr->cbegin(), args.ICs.size() - 1); 
      auto end = std::prev(args.time_mesh_ptr->cend());

      exec.ConsumeSolutionList(args.ICs.begin(), args.ICs.end()); // no more .cbegin() -> should use move semantics
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), std::next(it)); 

      double t = *it; 
      while(it!= end)
      {
        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ 
            ((lam_args.template BeforeLinAlgebra<decltype(it),decltype(exec),decltype(m_rhs),OSteps::FDStep_Type::IMPLICIT>(it, args.domain_mesh_ptr, exec, m_rhs)), ...); 
          }, 
          m_ostep_tup
        ); 

        Mat = m_rhs.GetMat(); 

        // working on right end of [t(n-1), t(n)]
        ++it;
        t = *it;  

        exec.BuildNextTime(t); 

        // scale Mat according to 1 / dt ... 
        if constexpr(std::is_same<decltype(exec.inv_coeff()), const double&>::value){
          // INV_COEFF_T is a scalar
          Mat  *= exec.inv_coeff(); 
        }
        else{
          // INV_COEFF_T is a Matrix 
          TExprs::MatrixStorage_t temp = exec.inv_coeff() * Mat; 
          Mat = std::move(temp); 
        }

        // outside steps matrix before step(Mat) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, Mat)), ...); }, 
          m_ostep_tup
        ); 

        Eigen::VectorXd rhs = std::move(exec.RhsVector()); 
        // outside steps solution before step (rhs) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolBeforeStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, rhs)), ...); }, 
          m_ostep_tup
        ); 

        // Explicit Step 
        Eigen::VectorXd next_sol = Mat * exec.MostRecentSol() + rhs; 
        
        // ourside steps solution after step(next_sol) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolAfterStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, next_sol)), ...); }, 
          m_ostep_tup
        ); 

        // give expiring solution to WRITE_POLICY_T
        write_policy.SaveSolution(std::move(exec.ExpiringSol())); 

        // push Solution, time to executor 
        exec.ConsumeSolution(next_sol);
        exec.ConsumeTime(t); 
      }    

      // Write remaining solutions to write_policy 
      for(auto i=0; i<exec.StoredSols().size()-1; i++) write_policy.SaveSolution( std::move(exec.StoredSols()[i])); 

      // write policy also determines return type / how to handle last solution
      return write_policy.ConsumeLastSolution(std::move( exec.MostRecentSol() )); 
    
    } // end .Calculate(args, write_policy) 

}; 

} // end namespace Solvers 

#endif // ExplicitSolver.hpp 