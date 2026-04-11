// OStepBase.hpp 
//
// CRTP bases for 1D / XD Outside Step. 
// Boundary Conditions will be viewed as another type of OStep in the new framework 
//
// JAF 1/31/2026 

#ifndef OSTEPBASE_H
#define OSTEPBASE_H

#include<LinOps/LinearOpBase.hpp> // LinOps::MatrixStorage_t + forward declarations 

namespace OSteps{
using StridedRef_t = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
using LinOps::MatrixStorage_t; 
using LinOps::SharedConstMesh1D; 
using LinOps::SharedConstMeshXD;

enum class FDStep_Type{
  EXPLICIT, 
  IMPLICIT 
}; 

enum class StepDims_Flag{
  OneDim, 
  XDim 
}; 

template<typename DERIVED, StepDims_Flag DIMS_FLAG>
class OStepBase
{
  public:
    // Type Defs + Compile time flags -----------------
    constexpr static StepDims_Flag dims_flag = DIMS_FLAG;
    // constexpr static bool is_time_dep_flag = false; // going to remove in favor of SetTime outside step ... 
    using MESH_T = std::conditional_t<DIMS_FLAG == StepDims_Flag::OneDim, SharedConstMesh1D, SharedConstMeshXD>; 

    // Member Funcs ------------------------------------------------
    template<typename TIME_ITER, typename LHS_EXECUTOR, typename RHS_EXPR, FDStep_Type Step>
    void BeforeLinAlgebra(TIME_ITER& time_iter, MESH_T& mesh, LHS_EXECUTOR& exec, RHS_EXPR& rhs_expr)
    {/* this is defaulted to do nothing.*/}

    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void MatBeforeStep(double t, const MESH_T& mesh, LinOps::MatrixStorage_t& Mat) const 
    {
      static_cast<const DERIVED*>(this)->template MatBeforeStep<STEP>(t,mesh, Mat); 
    }
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void SolBeforeStep(double t, const MESH_T& mesh, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->template SolBeforeStep<STEP>(t,mesh, Sol); 
    }
    template<FDStep_Type STEP = FDStep_Type::EXPLICIT>
    void SolAfterStep(double t, const MESH_T& mesh, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->template SolAfterStep<STEP>(t,mesh, Sol); 
    }
}; 

template<typename DERIVED>
using OStepBase1D = OStepBase<DERIVED, StepDims_Flag::OneDim>;

template<typename DERIVED>
using OStepBaseXD = OStepBase<DERIVED, StepDims_Flag::XDim>;

} // end namesace OSteps

#endif // OStepBase.hpp