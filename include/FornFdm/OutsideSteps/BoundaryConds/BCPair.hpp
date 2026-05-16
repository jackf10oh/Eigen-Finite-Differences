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
  template<typename... BCPairs_Ts>
  friend class BCList; 
  private:
    // Member Data -----------------------------------------------------------
    typename std::remove_reference<LBC_T>::type m_left; 
    typename std::remove_reference<RBC_T>::type m_right; 
    
  public:
    // Constructors + Destructor =================================================
    BCPair() = delete;

    BCPair(LBC_T l, RBC_T r)
      : m_left(l),m_right(r)
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
        m_left.SetStencilL(t.next, m->getAxis(0), Mat); 
        m_right.SetStencilR(t.next, m->getAxis(0), Mat); 
      }
    }

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fornfdm::StridedRef u, const TCtx& t, const Ctx& ctx) const
    {
      auto m = ctx.getMesh(); 
      if (m->numDims() != 1) throw std::runtime_error("incorrect # of dims passed to 1D boundary condition");        
      if constexpr(STEP == StepType::Implicit){
        m_left.SetImpSolL(t.next, m->getAxis(0), u); 
        m_right.SetImpSolR(t.next, m->getAxis(0), u); 
      }
    }

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fornfdm::StridedRef u, const TCtx& t, const Ctx& ctx) const 
    {
      auto m = ctx.getMesh(); 
      if (m->numDims() != 1) throw std::runtime_error("incorrect # of dims passed to 1D boundary condition"); 
      if constexpr(STEP == StepType::Explicit){
        m_left.SetSolL(t.next, m->getAxis(0), u); 
        m_right.SetSolR(t.next, m->getAxis(0), u); 
      }  
    }
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // BCPair.hpp