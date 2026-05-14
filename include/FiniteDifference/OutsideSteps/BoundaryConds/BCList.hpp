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

template<typename... BCPairs_Ts>
class BCList : public OStepBase<BCList<BCPairs_Ts...>>
{
  private:
    // Type Traits + Utils ----------------
    template<typename T>
    struct is_bc_pair_impl : public std::false_type{}; 

    template<typename L, typename R>
    struct is_bc_pair_impl<BCPair<L,R>> : public std::true_type{}; 

    template<typename T>
    using is_bc_pair = is_bc_pair_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

    template<typename T, std::size_t ignore=0>
    using repeat_type = T; 

    template<typename T, std::size_t... Is>
    static auto make_repeat_tuple(std::index_sequence<Is...>)
    {
      return std::tuple< repeat_type<T,Is>... >{}; 
    }
    // member data. -------------------------------------------------
    // list of boundary conditions. 1 per Dimension 
    std::tuple< std::remove_reference_t<BCPairs_Ts>... > m_bcs_list; 
    using MATS_T = decltype(make_repeat_tuple<fornfdm::CSRMatrix>(std::make_index_sequence<sizeof...(BCPairs_Ts)>{})); 
    MATS_T m_mats; 
    
  public:
    // Constructors + Destructors =========================================
    BCList()=delete; 
    
    template<typename = std::enable_if_t<
      std::conjunction_v<
          is_bc_pair<BCPairs_Ts> ...
        >
      >
    >
    BCList(BCPairs_Ts... args) 
      : m_bcs_list(args...), 
      m_mats(make_repeat_tuple<fornfdm::CSRMatrix>(std::make_index_sequence<sizeof...(BCPairs_Ts)>{}))
    { 
      // reserves 10 entries inside of a fornfdm::Matrix
      auto reserve_lam = [](auto& mat){mat.reserve(10);}; 

      // go through each fornfdm::Matrix and reserve 10 entries
      std::apply(
        [&](auto&... mats){ ((reserve_lam(mats)), ...); }, 
        m_mats
      ); 
    } 
    
    BCList(const BCList& other)=default;
    
    // destructor 
    ~BCList()=default; 

    // Member Funcs =================================================
    template<StepType STEP, typename TIMECTX = TimeContext<>, typename CONSTCTX = Context<> >
    void MatBeforeStep(fornfdm::CSRMatrix& Mat, const TIMECTX& t = {}, const CONSTCTX& ctx = {})
    {      
      // check args are compaitble 
      assert((ctx.getMesh()->numDims() == sizeof...(BCPairs_Ts)) && "BCList MatBeforeStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == StepType::Implicit)
      {    
        constexpr std::size_t N = sizeof...(BCPairs_Ts); 
        // Sets m_mats 2nd to last according to bc pairs 2nd to last... 
        prepare_flat_stencils(t.next, ctx.getMesh(), std::make_index_sequence<N>{}); 

        // resize + set first matrix in m_mats 
        std::get<0>(m_mats).resize(ctx.getMesh()->sizeOfDim(0), ctx.getMesh()->sizeOfDim(0)); 
        std::get<0>(m_mats).setZero();
        std::get<0>(m_bcs_list).m_left.SetStencilL(t.next,ctx.getMesh()->getAxis(0), std::get<0>(m_mats));  
        std::get<0>(m_bcs_list).m_right.SetStencilR(t.next,ctx.getMesh()->getAxis(0), std::get<0>(m_mats));  

        decltype(auto) mask = make_overwrite_mask(ctx.getMesh(), std::make_index_sequence<N>{}); 
        fornfdm::utils::overwrite_stencil(Mat, mask); 
      }
    } // end MatBeforeStep 

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecBeforeStep(fornfdm::StridedRef u, const TCtx& t, const Ctx& ctx)
    {
      // check args are compaitble 
      auto mesh = ctx.getMesh(); 
      assert((mesh->numDims() == sizeof...(BCPairs_Ts)) && "BCList VecBeforeStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == StepType::Implicit)
      {
        // this lambda takes 1 BCPair<L,R> and applies it to Sol
        // without double assigning to corners/edges of XDim space 

        auto set_dim_boundaries_imp = [&](const auto& bc_pair, std::size_t dim){
          // Mesh1D that this bc_pair operates on  
          const auto& mesh_1dim = mesh->getAxis(dim); 
          // vector of eigen stride views that "look" like Discretization1Ds along mesh_1dim 
          using Stride_t = typename Eigen::Stride<0,Eigen::Dynamic>; 
          using StrideView_t = typename Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t>; 
          std::size_t mod = mesh->sizesMiddleProduct(0,dim); 
          std::size_t num_copies = mod * mesh->sizesMiddleProduct(dim+1, mesh->numDims()); 
          std::size_t s0 = mesh->sizeOfDim(dim); 
          std::size_t scale = mod * s0;  
          Stride_t stride(0, mod); 
          
          // iterate through views that look like Mesh1D
          for(std::size_t i=0; i < num_copies; i++){
            // determine if it has been set by lower dimension 
            bool set_by_low_dim = false; 
            std::size_t s1 = 1; 
            for(int dim_i=0; dim_i<dim; dim_i++){
              std::size_t s2 = mesh->sizeOfDim(dim_i); 
              // in first group of values
              if((i/s1)%s2 == 0) set_by_low_dim = true;  
              // in last group of values
              if((i/s1)%s2 == s2-1) set_by_low_dim = true;  
              // next check uses a large bucket 
              s1 *= s2; 
            }
            // if it hasn't, use BC on it. 
            if(!set_by_low_dim){
              std::size_t offset = (i % mod) + (scale * (i/mod)); // note: was previously using (mod? i%mod : i). shouldn't need it if mod is never == 0
              StrideView_t view(u.data()+offset, s0, stride); 
              bc_pair.m_left.SetImpSolL(t.next, mesh_1dim, view); 
              bc_pair.m_right.SetImpSolR(t.next, mesh_1dim, view); 
            }
          }
        }; // end set_dim_boundaries lambda 
        
        // apply set_dim_boundaries lambda along each dimension 
        std::size_t dim = 0; 
        std::apply(
          [&](const auto&... args){
            (set_dim_boundaries_imp(args, dim++), ...);
          }, 
          m_bcs_list
        ); 
      } // end if constexpr(step)
    } // end VecBeforeStep 

    template<StepType STEP, typename TCtx=TimeContext<>, typename Ctx=Context<> >
    void VecAfterStep(fornfdm::StridedRef u, const TCtx& t, const Ctx& ctx)
    {
      
      auto mesh = ctx.getMesh(); // get the MeshXD 
      assert((mesh->numDims() == sizeof...(BCPairs_Ts)) && "BCList VecAfterStep error: Mesh.numDims() != size of tuple / list of 1D BCs "); 
      
      if constexpr(STEP == StepType::Explicit)
      {
        // this lambda takes 1 BCPair<L,R> and applies it to Sol
        // without double assigning to corners/edges of XDim space 

        auto set_dim_boundaries = [&](const auto& bc_pair, std::size_t dim){
          using Stride_t = typename Eigen::Stride<0,Eigen::Dynamic>; 
          using StrideView_t = typename Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t>; 
          // Mesh1D that this bc_pair operates on  
          const auto& mesh_1dim = mesh->getAxis(dim); // get the Mesh1D 
          // vector of eigen stride views that "look" like Discretization1Ds along mesh_1dim 
          std::size_t mod = mesh->sizesMiddleProduct(0,dim); 
          std::size_t num_copies = mod * mesh->sizesMiddleProduct(dim+1, mesh->numDims()); 
          std::size_t s0 = mesh->sizeOfDim(dim); 
          std::size_t scale = mod * s0;  
          Stride_t stride(0, mod); 

          // iterate through views that look like Mesh1D
          for(std::size_t i=0; i < num_copies; i++){
            // determine if it has been set by lower dimension 
            bool set_by_low_dim = false; 
            std::size_t s1 = 1; 
            for(int dim_i=0; dim_i<dim; dim_i++){
              std::size_t s2 = mesh->sizeOfDim(dim_i); 
              // in first group of values
              if((i/s1)%s2 == 0) set_by_low_dim = true;  
              // in last group of values
              if((i/s1)%s2 == s2-1) set_by_low_dim = true;  
              // next check uses a large bucket 
              s1 *= s2; 
            }
            // if it hasn't, use BC on it. 
            if(!set_by_low_dim){
              std::size_t offset = (mod ? i % mod : i) + (scale * (i/mod));
              StrideView_t view(u.data()+offset, s0, stride); 
              bc_pair.m_left.SetSolL(t.next, mesh_1dim, view); 
              bc_pair.m_right.SetSolR(t.next, mesh_1dim, view); 
            }
          } // end for loop through 1 dimensional views 

        }; // end set_dim_boundaries lambda 

        // apply set_dim_boundaries lambda along each dimension 
        std::size_t dim = 0; 
        std::apply(
          [&](const auto&... args){
            (set_dim_boundaries(args, dim++), ...);
          }, 
          m_bcs_list
        ); 
      } // end if constexpr(step)
    } // end VecAfterStep

  private:
    // Unreachable =========================================================== 
    template<std::size_t Idx>
    void flat_stencil(fornfdm::Real t, const fornfdm::Mesh* m)
    {
      // !!! uses first entry in m_mats as temp storage !!! 
      if constexpr( Idx > 0 ){
        // using ith 1D mesh out of MeshXD 
        const auto& mesh = m->getAxis(Idx); 

        // resize 1st matrix + ith matrix to 1 x N where N == mesh.size() 
        std::get<0>(m_mats).resize(1,mesh.size()); 
        std::get<0>(m_mats).setZero(); 
        std::get<Idx>(m_mats).resize(1,mesh.size()); 
        std::get<Idx>(m_mats).setZero(); 

        // set ith matrix to be top row of stencil, 1st matrix to bottom row of stencil 
        std::get<Idx>(m_bcs_list).m_left.SetStencilL(t,mesh,std::get<Idx>(m_mats)); 
        std::get<Idx>(m_bcs_list).m_right.SetStencilR(t,mesh,std::get<0>(m_mats));

        // copy values from "bottom" row into "top" row. 
        fornfdm::CSRMatrix::InnerIterator it(std::get<0>(m_mats),0); 
        for(; it; ++it) std::get<Idx>(m_mats).coeffRef(it.row(),it.col()) = it.value();  
      }
    }

    template<std::size_t... Is>
    void prepare_flat_stencils(fornfdm::Real t, const fornfdm::Mesh* m, std::index_sequence<Is...>){
      (flat_stencil<Is>(t,m), ...); 
    } 

    template<std::size_t Idx>
    auto highdim_stencil(const fornfdm::Mesh* m)
    {
      if constexpr(Idx == 0){
        // first matrix just gets kronecker product into high dim.
        return fornfdm::utils::make_BlockDiag(std::get<Idx>(m_mats), m->sizesMiddleProduct(1, sizeof...(BCPairs_Ts))); 
      }
      else{
        // other dimensions have to be taken 
        // flat stencil -> sparse diag (repeats) -> kronecker product 
        std::size_t s1 = m->sizesMiddleProduct(0, Idx);
        std::size_t s2 = m->sizesMiddleProduct(Idx+1, sizeof...(BCPairs_Ts)); 
        return fornfdm::utils::make_BlockDiag(fornfdm::utils::make_SparseDiag<fornfdm::utils::SparseDiagPattern::REPEAT>(std::get<Idx>(m_mats), s1), s2); 
      }
    }

    template<std::size_t... Is>
    decltype(auto) make_overwrite_mask(const fornfdm::Mesh* m, std::index_sequence<Is...>)
    {
      if constexpr(sizeof...(Is)==1)
      {
        return std::as_const(std::get<0>(m_mats)); 
      }
      else
      {
        return fornfdm::utils::make_FillRow_fold( highdim_stencil<Is>(m) ...); 
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