// OStepBase.hpp 
//
// CRTP bases for Outside Step. 
// Boundary Conditions will be viewed as another type of OStep in the new framework 
//
// JAF 1/31/2026 

#ifndef FORNFDM_OSTEPS_OSTEPBASE_H
#define FORNFDM_OSTEPS_OSTEPBASE_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include "../types.hpp"
#include "StepContexts.hpp"

namespace fornfdm{
  namespace osteps{    

enum class StepType{
  Explicit, 
  Implicit 
}; 

template<typename DERIVED>
struct OStepBase
{
  public:
    // Member Funcs ------------------------------------------------
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void beforeLinAlgebra(const TCtx& t, Ctx& ctx){/* edit everything before linear algebra starts*/}

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void applyBeforeMat(fornfdm::CSRMatrix& Mat, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}
    
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void applyBeforeVec(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void applyAfterVec(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx){/* edit the solution vector after the step */}

}; 

  } // end namespace osteps 
} // end namespace fornfdm 

#endif // OStepBase.hpp