// SetTimeOstep.hpp
//
// Outside step that only has an effect on LHS TExpr executor + RHS LinOp expression 
// by calling their SetTime() member functions 
// 
// JAF 2/3/2026 

#ifndef SETTIMEOSTEP_H
#define SETTIMEOSTEP_H

#include "OStepBase.hpp"
#include "StepContexts.hpp"

namespace OSteps{

class SetTimeOStep: OStepBase<SteTimeOStep>
{
  public:
    // Member Funcs ------------------------------------------------
    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void BeforeLinAlgebra(const TCtx& t, Ctx& ctx)
    {
      // Calls SetTime on any LinOps inside the expressions 
      if constexpr(Step == FDStep_Type::EXPLICIT){
        ctx.getExecutor().expr_SetTime(t.now); 
        ctx.getRhsExpr().SetTime(t.now); 
      }
      if constexpr(Step == FDStep_Type::IMPLICIT){
        ctx.getExecutor().expr_SetTime(t.next); 
        ctx.getRhsExpr().SetTime(t.next);
      }
    }
}; 

} // end namespace OSteps

#endif 