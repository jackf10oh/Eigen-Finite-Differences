// OStepBase.hpp 
//
// CRTP bases for 1D / XD Outside Step. 
// Boundary Conditions will be viewed as another type of OStep in the new framework 
//
// JAF 1/31/2026 

#ifndef OSTEPBASE_H
#define OSTEPBASE_H

#include "../LinOps/LinOpTraits.hpp" 

#include "StepContexts.hpp"

namespace fdm{
  namespace osteps{    

enum class StepType{
  Explicit, 
  Implicit 
}; 

template<typename DERIVED>
class OStepBase
{
  public:
    // Member Funcs ------------------------------------------------
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void BeforeLinAlgebra(const TCtx& t, Ctx& ctx){/* edit everything before linear algebra starts*/}

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void MatBeforeStep(fdm::Matrix& Mat, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}
    
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fdm::StridedRef u, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fdm::StridedRef u, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}

}; 

  } // end namespace osteps 
} // end namespace fdm 

#endif // OStepBase.hpp