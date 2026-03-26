// OStepBase.hpp 
//
// CRTP bases for 1D / XD Outside Step. 
// Boundary Conditions will be viewed as another type of OStep in the new framework 
//
// JAF 1/31/2026 

#ifndef OSTEPBASE_H
#define OSTEPBASE_H

#include<LinOps/LinearOpBase.hpp> // LinOps::MatrixStorage_t + forward declarations 

#include "StepContexts.hpp"

namespace OSteps{

using StridedRef_t = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
using LinOps::MatrixStorage_t;

enum class FDStep_Type{
  EXPLICIT, 
  IMPLICIT 
}; 

template<typename DERIVED>
class OStepBase
{
  public:
    // Member Funcs ------------------------------------------------
    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void BeforeLinAlgebra(const TCtx& t, Ctx& ctx){/* edit everything before linear algebra starts*/}

    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void MatBeforeStep(LinOps::MatrixStorage_t& Mat, const TCtx& t, const Ctx& ctx) const {/* edit the solution vector after the step */}
    
    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(StridedRef_t u, const TCtx& t, const Ctx& ctx) const {/* edit the solution vector after the step */}

    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(StridedRef_t u, const TCtx& t, const Ctx& ctx) const {/* edit the solution vector after the step */}

}; 

} // end namesace OSteps

#endif // OStepBase.hpp