// ForcingTerm.hpp
//
// represents f(t,x,y,z) in equations like 
// Ut = Ux + f(t,x,y,z) 
// will always be added to vector before step.
// i.e. results in implicit explicit (IMEX) schemes.   
// typically placed before boundary conditions in osteps tuples
//
// JAF 5/7/2026 

#ifndef FORNFDM_OSTEPS_FORCINGTERM_H
#define FORNFDM_OSTEPS_FORCINGTERM_H

#include "../traits.hpp" // callable_traits 
#include "../Coordinate.hpp"
#include "../Mesh.hpp" 
#include "OStepBase.hpp" 
#include "StepContexts.hpp"

namespace fornfdm{
namespace osteps{

template<class Callable>
class ForcingTerm : public OStepBase<ForcingTerm<Callable>>
{
  private:
    // Type Defs -------------------- 
    using CallableCleaned = std::remove_reference_t<Callable>; 
    // Member Data --------------------- 
    CallableCleaned m_functor; 

  public:
    // Constructors + Destructor ====================
    ForcingTerm()=delete; 

    ForcingTerm(Callable f)
      : m_functor(f)
    {}

    ForcingTerm(const ForcingTerm& other)=default; 

    // destructor 
    ~ForcingTerm()=default; 

    // Member Functions ---------------- 
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void applyBeforeVec(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx) const 
    {
      using Binded = typename fornfdm::internal::BindFirst<CallableCleaned>; 
      if constexpr(STEP==StepType::Explicit)
      {
        // use left end point t[n]
        auto xpr = fornfdm::discretize(ctx.getMesh(), Binded(m_functor, t.now)); 
        u += ctx.getExecutor()->getInvCoeff() * xpr; 
      }
      else if constexpr(STEP==StepType::Implicit)
      {
        // use midpoint (t[n] + t[n+1]) / 2  
        auto xpr = fornfdm::discretize(ctx.getMesh(), Binded(m_functor, (t.now+t.next)/2));
        u += ctx.getExecutor()->getInvCoeff() * xpr; 
      }
    }

    // getter utility to stored callable. 
    const auto& functor() const { return m_functor; }

}; 

} // end namespace osteps 
} // end namespace fornfdm 

#endif // ForcingTerm.hpp 