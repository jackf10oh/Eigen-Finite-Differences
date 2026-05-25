// BCPair.hpp
//
// Holds 2 BCs that are applied 
// to left/right side of vector 
// and top/bottom rows of CSRMatrix 
//
// JAF 12/8/2025

#ifndef FORNFDM_OSTEPS_BCPAIR_H
#define FORNFDM_OSTEPS_BCPAIR_H

#include "../../Types.hpp"
#include "../OStepBase.hpp"

namespace fornfdm{
  namespace osteps{

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
template<typename LBC_T,typename RBC_T>
class BCPair: public OStepBase<BCPair<LBC_T,RBC_T>>
{
  public:
    // Member Data -----------------------------------------------------------
    typename std::remove_reference<LBC_T>::type left_bc; 
    typename std::remove_reference<RBC_T>::type right_bc; 
    
    // Constructors + Destructor =================================================
    BCPair() = delete;

    BCPair(LBC_T l, RBC_T r)
      : left_bc(l),right_bc(r)
    {};

    BCPair(const BCPair& other)=default; 

    // destructor
    virtual ~BCPair()=default; 

    // Member Functions ==================================================================
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void MatBeforeStep(fornfdm::CSRMatrix& Mat, const TCtx& t, const Ctx& ctx) const
    {
      auto m = ctx.getMesh(); 
      if (m->numDims() != 1) throw std::runtime_error("incorrect # of dims passed to 1D boundary condition"); 
      if constexpr(STEP == StepType::Implicit){ 
        left_bc.SetStencilL(t.next, m->getAxis(0), Mat); 
        right_bc.SetStencilR(t.next, m->getAxis(0), Mat); 
      }
    }

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx) const
    {
      auto m = ctx.getMesh(); 
      if (m->numDims() != 1) throw std::runtime_error("incorrect # of dims passed to 1D boundary condition");        
      if constexpr(STEP == StepType::Implicit){
        left_bc.SetImpSolL(t.next, m->getAxis(0), u); 
        right_bc.SetImpSolR(t.next, m->getAxis(0), u); 
      }
    }

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx) const 
    {
      auto m = ctx.getMesh(); 
      if (m->numDims() != 1) throw std::runtime_error("incorrect # of dims passed to 1D boundary condition"); 
      if constexpr(STEP == StepType::Explicit){
        left_bc.SetSolL(t.next, m->getAxis(0), u); 
        right_bc.SetSolR(t.next, m->getAxis(0), u); 
      }  
    }
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // BCPair.hpp