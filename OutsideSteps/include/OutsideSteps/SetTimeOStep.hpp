// SetTimeOstep.hpp
//
// Outside step that only has an effect on LHS TExpr executor + RHS LinOp expression 
// by calling their setTime() member functions 
// 
// JAF 2/3/2026 

#ifndef SETTIMEOSTEP_H
#define SETTIMEOSTEP_H

#include "OStepBase.hpp"

namespace OSteps{

class SetTimeOStep: public OStepBaseXD<>, public OStepBase1D<>
{
  public:
    // Calls setTime on any LinOps inside the expression
    template<typename TIME_ITER, typename LHS_EXECUTOR, typename RHS_EXPR, FDStep_Type Step>
    void BeforeLinAlgebra(TIME_ITER& time_iter, Mesh1D_SPtr_t& mesh, LHS_EXECUTOR& exec, RHS_EXPR& rhs_expr) 
    {
      if constexpr(Step == FDStep_Type::EXPLICIT){
        exec.expr_SetTime(*time_iter); 
        rhs_expr.setTime(*time_iter); 
        exec.set_mesh(mesh); 
        rhs_expr.set_mesh(mesh); 
      }
      if constexpr(Step == FDStep_Type::IMPLICIT){
        exec.expr_SetTime(*(time_iter+1)); 
        rhs_expr.setTime(*(time_iter+1)); 
        exec.set_mesh(mesh); 
        rhs_expr.set_mesh(mesh); 
      }
    }

    template<typename TIME_ITER, typename LHS_EXECUTOR, typename RHS_EXPR, FDStep_Type Step>
    void BeforeLinAlgebra(TIME_ITER& time_iter, MeshXD_SPtr_t& mesh, LHS_EXECUTOR& exec, RHS_EXPR& rhs_expr) 
    {
      if constexpr(Step == FDStep_Type::EXPLICIT){
        exec.expr_SetTime(*time_iter); 
        rhs_expr.setTime(*time_iter); 
        exec.set_mesh(mesh); 
        rhs_expr.set_mesh(mesh); 
      }
      if constexpr(Step == FDStep_Type::IMPLICIT){
        exec.expr_SetTime(*(time_iter+1)); 
        rhs_expr.setTime(*(time_iter+1)); 
        exec.set_mesh(mesh); 
        rhs_expr.set_mesh(mesh); 
      }
    }


    // all of the other member functions don't change anything 
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void MatBeforeStep(double t, const LinOps::Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat) const {}
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void SolBeforeStep(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const {}
    template<FDStep_Type STEP = FDStep_Type::EXPLICIT>
    void SolAfterStep(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const {} 
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void MatBeforeStep(double t, const LinOps::MeshXD_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat) const {}
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void SolBeforeStep(double t, const LinOps::MeshXD_SPtr_t& mesh, StridedRef_t Sol) const {}
    template<FDStep_Type STEP = FDStep_Type::EXPLICIT>
    void SolAfterStep(double t, const LinOps::MeshXD_SPtr_t& mesh, StridedRef_t Sol) const {} 
}

} // end namespace OSteps

#endif 