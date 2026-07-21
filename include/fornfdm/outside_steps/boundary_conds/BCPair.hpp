// BCPair.hpp
//
// Holds 2 BCs that are applied 
// to left/right side of vector 
// and top/bottom rows of CSRMatrix 
//
// JAF 12/8/2025

#ifndef FORNFDM_OSTEPS_BCPAIR_H
#define FORNFDM_OSTEPS_BCPAIR_H

#include<array>
#include<cstdint>
#include<type_traits>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include "../../types.hpp"
#include "../OStepBase.hpp"

namespace fornfdm{
  namespace osteps{

// Utility class for overwriting top/bottom rows in a stencil
struct BoundaryRow
{
  std::array<fornfdm::Scalar,3> vals; 
  std::size_t nnz;
};

// Base Class for Boundary Conditions. all operators make no changes to stencil / solution 
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
    ~BCPair()=default; 

    // Member Functions ==================================================================
    using OStepBase<BCPair<LBC_T,RBC_T>>::applyBeforeMat;
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    std::enable_if_t<STEP==StepType::Implicit, void> applyBeforeMat(fornfdm::CSRMatrix& Mat, const TCtx& t, const Ctx& ctx) const
    {
      auto m = ctx.getMesh(); 
      assert((m->numDims() == 1) && "incorrect # of dims passed to 1D boundary condition. use BCList instead");
      const auto& axis = m->getAxis(0);
      Mat.row(0) *= 0.0;
      BoundaryRow top_row = left_bc.getTopRow(t.next, axis);
      for(auto n=0; n<top_row.nnz; ++n)
      {
        Mat.coeffRef(0,n) = top_row.vals[n];
      }
      Mat.row(axis.size()-1) *= 0.0;
      BoundaryRow bottom_row = right_bc.getBottomRow(t.next, axis);
      for(auto n=0; n<bottom_row.nnz; ++n)
      {
        Mat.coeffRef(axis.size()-1,axis.size()-bottom_row.nnz+n) = bottom_row.vals[n];
      }
    }

    using OStepBase<BCPair<LBC_T,RBC_T>>::applyBeforeVec;
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    std::enable_if_t<STEP==StepType::Implicit, void> applyBeforeVec(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx) const
    {
      auto m = ctx.getMesh(); 
      assert((m->numDims()==1) && "incorrect # of dims passed to 1D boundary condition");        
      left_bc.setImpSolLeft(t.next, m->getAxis(0), u); 
      right_bc.setImpSolRight(t.next, m->getAxis(0), u); 
    }

    using OStepBase<BCPair<LBC_T,RBC_T>>::applyAfterVec;
    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    std::enable_if_t<STEP==StepType::Explicit, void> applyAfterVec(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx) const 
    {
      auto m = ctx.getMesh(); 
      assert((m->numDims() == 1) && "incorrect # of dims passed to 1D boundary condition"); 
      left_bc.setExpSolLeft(t.next, m->getAxis(0), u); 
      right_bc.setExpSolRight(t.next, m->getAxis(0), u); 
    }
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // BCPair.hpp