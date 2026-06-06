// BCList.hpp
//
// holds a tuple of BCLists and applies them 
// in order to each dimension 
//
// JAF 2/1/2026  

#ifndef FORNFDM_OSTEPS_BCLIST_H
#define FORNFDM_OSTEPS_BCLIST_H

#include<cassert>
#include<tuple>
#include<Eigen/Dense>
#include<utility>
#include "../../Mesh.hpp"
#include "../../Utilities/FillRowExpr.hpp"
#include "../../Utilities/FillStencil.hpp"
#include "../../Utilities/BlockDiagExpr.hpp"
#include "../../Utilities/SparseDiagExpr.hpp"
#include "../OStepBase.hpp"
#include "BCPair.hpp"

namespace fornfdm{
  namespace osteps{

template<typename... BCPairTypes>
class BCList : public OStepBase<BCList<BCPairTypes...>>
{
  public:
    // member data. -------------------------------------------------
    constexpr static std::size_t numDims = sizeof...(BCPairTypes);
    // list of boundary conditions. 1 per Dimension 
    std::tuple< std::remove_reference_t<BCPairTypes>... > bcs_list; 

  public:
    // Constructors + Destructors =========================================
    BCList()=delete; 
    
    BCList(BCPairTypes... args) 
      : bcs_list(args...)
    {} 
    
    BCList(const BCList& other)=default;
    
    // destructor 
    ~BCList()=default; 

    // Member Funcs =================================================
    template<StepType STEP, typename TIMECTX = TimeContext<>, typename CONSTCTX = Context<> >
    void MatBeforeStep(fornfdm::CSRMatrix& mat, const TIMECTX& t, const CONSTCTX& ctx)
    {      
      // check args are compaitble 
      assert((ctx.getMesh()->numDims() == this->numDims) && "BCList MatBeforeStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == StepType::Implicit)
      {    
        auto rolling_prods = make_rolling_prods(ctx.getMesh());
        auto boundary_row_pairs = make_row_pairs(t.next, ctx.getMesh());
        setStencilNoFill<numDims-1>(mat, 0, rolling_prods, boundary_row_pairs);
      }
    } // end MatBeforeStep 

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx)
    {
      // check args are compaitble 
      auto mesh = ctx.getMesh(); 
      assert((mesh->numDims() == this->numDims) && "BCList VecBeforeStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == StepType::Implicit)
      {
        std::array<std::size_t, numDims> rolling_prods = make_rolling_prods((ctx.getMesh())); 
        setSolNoFill<true, numDims-1>(u.data(), 0, rolling_prods, ctx.getMesh(), t.next);
      } // end if constexpr(step)
    } // end VecBeforeStep 

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fornfdm::StrideRef u, const TCtx& t, const Ctx& ctx)
    {
      auto mesh = ctx.getMesh(); // get the MeshXD 
      assert((mesh->numDims() == this->numDims) && "BCList VecAfterStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 
      
      if constexpr(STEP == StepType::Explicit)
      {
        std::array<std::size_t, numDims> rolling_prods = make_rolling_prods((ctx.getMesh())); 
        setSolNoFill<false, numDims-1>(u.data(), 0, rolling_prods, ctx.getMesh(), t.next);
      } // end if constexpr(step)
    } // end VecAfterStep

  private:
    // Unreachable =========================================================== 
    auto make_rolling_prods(const fornfdm::Mesh* m)
    {
      constexpr std::size_t N = sizeof...(BCPairTypes); 
      std::array<std::size_t, N> result;
      auto rolling = m->sizeOfDim(0);
      result[0] = rolling; 
      for(auto i=1; i < N; ++i)
      {
        rolling *= m->sizeOfDim(i);
        result[i] = rolling; 
      }
      return result;
    }

    template<std::size_t Idx>
    inline auto make_row_pair_single(fornfdm::Real t, const fornfdm::Vector& axis) const
    {
      return std::make_pair(
        std::get<Idx>(bcs_list).left_bc.getTopRow(t,axis), 
        std::get<Idx>(bcs_list).right_bc.getTopRow(t,axis)
      ); 
    }

    template<std::size_t... Idxs>
    auto make_row_pairs_impl(fornfdm::Real t, const fornfdm::Mesh* m, std::index_sequence<Idxs...>) const
    {
      std::array<std::pair<osteps::BoundaryRow,osteps::BoundaryRow>, sizeof...(Idxs)> result;
      (
        (result[Idxs] = make_row_pair_single<Idxs>(t,m->getAxis(Idxs))),
         ...
      );
      return result;
    }

    inline auto make_row_pairs(fornfdm::Real t, const fornfdm::Mesh* m) const
    {
      return make_row_pairs_impl(t, m, std::make_index_sequence<numDims>{});
    }

    template<std::size_t ithDim>
    void setStencilNoFill(fornfdm::CSRMatrix& mat, 
                          std::size_t offset,
                           const std::array<std::size_t,numDims>& rolling_prods, 
                           const std::array<std::pair<BoundaryRow,BoundaryRow>,numDims>& boundary_row_pairs) const
    {
      if constexpr(ithDim == 0)
      {
        // use a row iterator to fill the top rows. 
        mat.row(offset) *= 0.0;
        const auto& top_row = boundary_row_pairs[0].first;
        for(auto n=0; n < top_row.nnz; ++n)
        {
          mat.coeffRef(offset, offset+n) = top_row.vals[n];
        } 
        // use a row iterator to fill the bottom rows.
        mat.row(offset + rolling_prods[0]-1) *= 0.0;
        const auto& bottom_row = boundary_row_pairs[0].second;
        for(auto n=0; n < bottom_row.nnz; ++n)
        {
          mat.coeffRef(offset + rolling_prods[0]-1, offset + rolling_prods[0] - bottom_row.nnz + n) = bottom_row.vals[n];
        } 
        // return/goto parent dimension
      }
      else
      {
        // setStencilFill() using THIS dimensions row iterator
        setStencilFill<ithDim-1, ithDim, true>(mat, offset, rolling_prods, boundary_row_pairs);
        // fill middle blocks with setStencilNoFill() and a certain offset 
        for(auto m=1; m < rolling_prods[ithDim-1]-1; ++m)
        {
          setStencilNoFill<ithDim-1>(mat, offset + m*rolling_prods[ithDim-1], rolling_prods, boundary_row_pairs);
        }
        // setStencilFill() using THIS dimensions row iterator 
        setStencilFill<ithDim-1,ithDim, false>(mat,offset + rolling_prods[ithDim] - rolling_prods[ithDim-1], rolling_prods, boundary_row_pairs); 
      }
    }

    template<std::size_t ithDim, std::size_t parentDim, bool useTopRow>
    void setStencilFill(fornfdm::CSRMatrix& mat, 
                        std::size_t offset,
                        const std::array<std::size_t,numDims>& rolling_prods, 
                        const std::array<std::pair<BoundaryRow,BoundaryRow>,numDims>& boundary_row_pairs) const
    {
      if constexpr(ithDim == 0)
      {
        // use a row iterator to fill the top row 
        mat.row(offset) *= 0.0;
        const auto& top_row = boundary_row_pairs[0].first;
        for(auto n=0; n<top_row.nnz; ++n)
        {
          mat.coeffRef(offset,offset + n) = top_row.vals[n];
        }

        // use parent row iterator to fill middle rows 
        if constexpr(useTopRow)
        {
          const auto& middle_row = boundary_row_pairs[parentDim].first;
          for(auto n=1; n<rolling_prods[0]-1; ++n)
          {
            mat.row(offset+n) *= 0.0;
            for(auto m=0; m<middle_row.nnz; ++m)
            {
              mat.coeffRef(offset+n, offset+ n + m * rolling_prods[parentDim-1]) = middle_row.vals[m];
            }
          }
        }
        else
        {
          const auto& middle_row = boundary_row_pairs[parentDim].second;
          for(auto n=1; n<rolling_prods[0]-1; ++n)
          {
            mat.row(offset+n) *= 0.0;
            for(auto m=0; m<middle_row.nnz; ++m)
            {
              mat.coeffRef(offset+n, offset + n - (middle_row.nnz-m-1) * rolling_prods[parentDim-1]) = middle_row.vals[m];
            }
          }
        }

        // use a row iterator to fill the bottom row 
        mat.row(offset + rolling_prods[0]-1) *= 0.0;
        const auto& bottom_row = boundary_row_pairs[0].second;
        for(auto n=0; n < bottom_row.nnz; ++n)
        {
          mat.coeffRef(offset + rolling_prods[0]-1, offset + rolling_prods[0] - bottom_row.nnz + n) = bottom_row.vals[n];
        } 

        // return/goto parent dim
      }
      else
      {
        // setStencilFill() using THIS dimensions row iterator
        setStencilFill<ithDim-1, ithDim, true>(mat, offset, rolling_prods, boundary_row_pairs);
        // fill middle blocks with setStencilFill() using parent row 
        for(auto m=1; m < rolling_prods[ithDim-1]-1; ++m)
        {
          setStencilFill<ithDim-1, parentDim, useTopRow>(mat, offset + m*rolling_prods[ithDim-1], rolling_prods, boundary_row_pairs);
        }
        // setStencilFill() using THIS dimensions row iterator 
        setStencilFill<ithDim-1,ithDim, false>(mat,offset + rolling_prods[ithDim] - rolling_prods[ithDim-1], rolling_prods, boundary_row_pairs); 
      }
    }

    template<bool useImplicit, std::size_t ithDim>
    void setSolNoFill(fornfdm::Scalar* data, 
                       std::size_t offset,
                       const std::array<std::size_t,numDims>& rolling_prods,
                       const fornfdm::Mesh* mesh,
                       fornfdm::Real t) const
    {
      if constexpr(ithDim == 0)
      {
        fornfdm::StrideView sol_1d(data + offset, mesh->sizeOfDim(0), fornfdm::Stride(0,1));
        if constexpr(useImplicit)
        {
          std::get<0>(bcs_list).left_bc.SetImpSolL(t, mesh->getAxis(0), sol_1d);
          std::get<0>(bcs_list).right_bc.SetImpSolR(t, mesh->getAxis(0), sol_1d);
        }
        else
        {
          std::get<0>(bcs_list).left_bc.SetSolL(t, mesh->getAxis(0), sol_1d);
          std::get<0>(bcs_list).right_bc.SetSolR(t, mesh->getAxis(0), sol_1d);
        }
      }
      else
      {
        // fill top block with THIS dimension as parent 
        setSolFill<useImplicit, ithDim-1, ithDim, true>(data, offset, rolling_prods, mesh, t);
        // middle blocks don't need filled
        for(auto n=1; n < mesh->sizeOfDim(ithDim)-1; ++n)
        {
          setSolNoFill<useImplicit, ithDim-1>(data, offset + n * rolling_prods[ithDim-1], rolling_prods, mesh, t);
        }
        // fill bottom block with THIS dimension as parent
        setSolFill<useImplicit, ithDim-1, ithDim, false>(data, offset + rolling_prods[ithDim] - rolling_prods[ithDim-1], rolling_prods, mesh, t);
      }
    }

    template<bool useImplicit, std::size_t ithDim, std::size_t parentDim, bool useLeftSide>
    void setSolFill(fornfdm::Scalar* data, 
                       std::size_t offset,
                       const std::array<std::size_t,numDims>& rolling_prods,
                       const fornfdm::Mesh* mesh,
                       fornfdm::Real t) const
    {
      if constexpr(ithDim == 0)
      {
        fornfdm::StrideView sol_1d(data + offset, mesh->sizeOfDim(0), fornfdm::Stride(0,1));
        if constexpr(useImplicit)
        {
          std::get<0>(bcs_list).left_bc.SetImpSolL(t, mesh->getAxis(0), sol_1d);
        }
        else
        {
          std::get<0>(bcs_list).left_bc.SetSolL(t, mesh->getAxis(0), sol_1d);
        }
        for(std::size_t n=1; n < mesh->sizeOfDim(0)-1; ++n)
        {
          if constexpr(useLeftSide)
          {
            // offset is just n steps forward from data+offset
            fornfdm::StrideView sol_1d(data + offset + n, mesh->sizeOfDim(parentDim), fornfdm::Stride(0,rolling_prods[parentDim]));
            if constexpr(useImplicit)
            {
              std::get<parentDim>(bcs_list).left_bc.SetImpSolL(t, mesh->getAxis(parentDim), sol_1d);
            }
            else
            {
              std::get<parentDim>(bcs_list).left_bc.SetSolL(t, mesh->getAxis(parentDim), sol_1d);
            }
          }
          else
          {
            // offset is n steps forward from data+offset
            // then backwards sizeOfDim(parentDim)-1 * stride
            // making end of data_offset ptr be the exact entry at data+offset+n
            fornfdm::Scalar* data_offset = data + offset + n - (mesh->sizeOfDim(parentDim) - 1) * rolling_prods[parentDim]; 
            fornfdm::StrideView sol_1d(data_offset, mesh->sizeOfDim(parentDim), fornfdm::Stride(0,rolling_prods[parentDim]));
            if constexpr(useImplicit)
            {
              std::get<parentDim>(bcs_list).right_bc.SetImpSolR(t, mesh->getAxis(parentDim), sol_1d);
            }
            else
            {
              std::get<parentDim>(bcs_list).right_bc.SetSolR(t, mesh->getAxis(parentDim), sol_1d);
            }
          }
        }
        if constexpr(useImplicit)
        {
          std::get<0>(bcs_list).right_bc.SetImpSolR(t, mesh->getAxis(0), sol_1d);
        }
        else
        {
          std::get<0>(bcs_list).right_bc.SetSolR(t, mesh->getAxis(0), sol_1d);
        }
      }
      else
      {
        // fill top block with THIS dimension as parent 
        setSolFill<useImplicit, ithDim-1, ithDim, true>(data, offset, rolling_prods, mesh, t);
        // fill middle blocks by parent dimensions
        for(auto n=1; n < mesh->sizeOfDim(ithDim)-1; ++n)
        {
          setSolFill<useImplicit, ithDim-1, parentDim, useLeftSide>(data, offset + n * rolling_prods[ithDim-1], rolling_prods, mesh, t);
        }
        // fill bottom block with THIS dimension as parent
        setSolFill<useImplicit, ithDim-1, ithDim, false>(data, offset + rolling_prods[ithDim] - rolling_prods[ithDim-1], rolling_prods, mesh, t);
      }
    }

}; // end class BCList 

  } // end namespace osteps
} // end namespace fornfdm 

#endif // BCList.hpp 

/* // order of priority for BCListXD.list to be applied 
corners / edges of XDim space are given priority to whichever BC would cover it first 
0 is unaffected by BCs. 
1 is first BC in .list,
2 is 2nd,
...  

// e.g. in 2D --------------------------
1 1 1 1 1 1 1
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
1 1 1 1 1 1 1


// e.g. in 3D -----------------------

// slice of bottom 
1 1 1 1 1 1 1
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
1 1 1 1 1 1 1

// slice directly off bottom. (3 no longer applies)
1 1 1 1 1 1 1
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
1 1 1 1 1 1 1

// so on and so forth 

*/