// BCPair.hpp
//
// 
//
// JAF 12/8/2025

#ifndef BCPAIR_H
#define BCPAIR_H

#include "../../LinOps/LinOpTraits.hpp" // StrideView, SharedConstMesh1D, etc 
#include "../../Mesh1D.hpp" 
#include "../OStepBase.hpp"

namespace fdm{
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
    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void MatBeforeStep(fdm::Matrix& Mat, const TCtx& t, const Ctx& ctx) const
    {
      using M = std::remove_cv_t<std::remove_reference_t<decltype(ctx.getMesh())>>; 
      static_assert(std::is_same_v<M,std::shared_ptr<const fdm::Mesh1D>>, "error in OSteps::BCPair. Domain is not Mesh1D"); 
      if constexpr(STEP == FDStep_Type::IMPLICIT){ 
        m_left.SetStencilL(t.next, ctx.getMesh(), Mat); 
        m_right.SetStencilR(t.next, ctx.getMesh(), Mat); 
      }
    }

    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fdm::StridedRef u, const TCtx& t, const Ctx& ctx) const
    {
      using M = std::remove_cv_t<std::remove_reference_t<decltype(ctx.getMesh())>>; 
      static_assert(std::is_same_v<M,std::shared_ptr<const fdm::Mesh1D>>, "error in OSteps::BCPair. Domain is not Mesh1D");       if constexpr(STEP == FDStep_Type::IMPLICIT){
        m_left.SetImpSolL(t.next, ctx.getMesh(), u); 
        m_right.SetImpSolR(t.next, ctx.getMesh(), u); 
      }
    }

    template<FDStep_Type STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fdm::StridedRef u, const TCtx& t, const Ctx& ctx) const 
    {
      using M = std::remove_cv_t<std::remove_reference_t<decltype(ctx.getMesh())>>; 
      static_assert(std::is_same_v<M,std::shared_ptr<const fdm::Mesh1D>>, "error in OSteps::BCPair. Domain is not Mesh1D"); 
      if constexpr(STEP == FDStep_Type::EXPLICIT){
        m_left.SetSolL(t.next, ctx.getMesh(), u); 
        m_right.SetSolR(t.next, ctx.getMesh(), u); 
      }  
    }
};

  } // end namespace osteps
} // end namespace fdm 

#endif // BCPair.hpp