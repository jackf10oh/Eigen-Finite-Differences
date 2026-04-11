// BCList.hpp
//
//
//
// JAF 2/1/2026  

#ifndef BCLIST_H
#define BCLIST_H 

#include<iostream>
#include<vector>
#include<tuple>
#include<unsupported/Eigen/KroneckerProduct>
#include<Utilities/FillStencil.hpp>

#include "../OStepBase.hpp"
#include "../BoundaryConds1D/BCPair.hpp"

namespace OSteps{

template<typename T>
struct is_bc_pair_impl : public std::false_type{}; 

template<typename L, typename R>
struct is_bc_pair_impl<BCPair<L,R>> : public std::true_type{}; 

template<typename T>
using is_bc_pair = is_bc_pair_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 


template<typename... BCPairs_Ts>
class BCList : public OStepBaseXD<BCList<BCPairs_Ts...>>
{
  public:
    // // compile time inspections ------------------------------ 
    // constexpr static bool is_time_dep_flag = false; // TODO! && fold of m_list. 
    // member data. -------------------------------------------------
    // list of boundary conditions. 1 per Dimension 
    std::tuple< std::remove_reference_t<BCPairs_Ts>... > m_list; 

    // Constructors + Destructors =========================================
    BCList()=delete; 
    template<typename = std::enable_if_t<
      std::conjunction_v<
          is_bc_pair<BCPairs_Ts> ...
        >
      >
    >
    BCList(BCPairs_Ts... args) : m_list(std::forward<BCPairs_Ts>(args)...){}; 
    BCList(const BCList& other)=default;
    // destructor 
    virtual ~BCList()=default; 

    // Member Funcs ------------------------------------------------
    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void MatBeforeStep(double t, const SharedConstMeshXD& mesh, LinOps::MatrixStorage_t& Mat) const 
    {      
      // check args are compaitble 
      if(mesh->numDims() != sizeof...(BCPairs_Ts)) throw std::runtime_error("BCList SolAfterStep error: MeshXD.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == FDStep_Type::IMPLICIT)
      {
        // initializations ---------------------------------------------
        // s1 = ith numDims size, s2 = cumulative product of all numDims < ith dim
        std::size_t s1=mesh->sizeOfDim(0); 
        std::size_t s2 = 1; 

        // lambda to give an NxN identity matrix 
        MatrixStorage_t cap_mat; 
        auto I = [&cap_mat](std::size_t N) -> const MatrixStorage_t& {cap_mat.resize(N,N); cap_mat.setIdentity(); return cap_mat;}; 
        // resulting mask for replacing rows in Mat input ---------------
        MatrixStorage_t mask;  
        
        // lambda to fill in rows of mask. leaving inner rows that aren't boundarys zero. 
        auto fill_mask_lambda = [&](const auto& bc_pair, std::size_t ith_dim)
        {
          if(ith_dim == 0)
          {
            mask = MatrixStorage_t(s1,s1); 
            bc_pair.template MatBeforeStep<FDStep_Type::IMPLICIT>(t,mesh->getMesh1D(ith_dim), mask); 
            s2 *= s1; 
            return; 
          }
          else
          {
            s1=mesh->sizeOfDim(ith_dim); 

            // take previous mask to higher dimension 
            MatrixStorage_t mask_temp = Eigen::KroneckerProductSparse(I(s1), mask); 
            
            // fill new empty rows with ith bc_pair
            MatrixStorage_t sp_diag = make_SparseDiag( flat_stencil(t, mesh->getMesh1D(ith_dim), bc_pair) , s2); 
            fill_stencil(mask_temp, sp_diag); 

            // take ownership of higher dim temp with mask 
            mask = std::move(mask_temp); 

            s2 *= s1; 
          }
        }; // end fill_mask_lambda 
        
        // loop through m_list tuple. 
        std::size_t ith_dim = 0; 
        std::apply(
          [&](const auto&... args){
            (fill_mask_lambda(args, ith_dim++), ...); 
          },
          m_list
        ); 
        // overwrite rows in Mat with nonzero rows from mask 
        overwrite_stencil(Mat, mask); 

      } // end if constexpr(step) 
    } // end MatBeforeStep<>(t,m,Mat) 

    template<FDStep_Type STEP = FDStep_Type::IMPLICIT>
    void SolBeforeStep(double t, const SharedConstMeshXD& mesh, StridedRef_t Sol) const 
    {
      // check args are compaitble 
      if(mesh->numDims() != sizeof...(BCPairs_Ts)) throw std::runtime_error("BCList SolAfterStep error: MeshXD.numDims() != size of tuple / list of 1D BCs "); 

      if constexpr(STEP == FDStep_Type::IMPLICIT)
      {
        // this lambda takes 1 BCPair<L,R> and applies it to Sol
        // without double assigning to corners/edges of XDim space 

        auto set_dim_boundaries_imp = [&](const auto& bc_pair, std::size_t dim){
          // Mesh1D that this bc_pair operates on  
          const auto& mesh_1dim = mesh->getMesh1D(dim); 
          // vector of eigen stride views that "look" like Discretization1Ds along mesh_1dim 
          auto views = mesh->makeOneDimViews(Sol, dim); 
          // iterate through views that look like Mesh1D
          for(std::size_t i=0; i<views.size(); i++){
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
              bc_pair.m_left.SetImpSolL(t, mesh_1dim, views[i]); 
              bc_pair.m_right.SetImpSolR(t, mesh_1dim, views[i]); 
            }
          }
        }; // end set_dim_boundaries lambda 
        
        // apply set_dim_boundaries lambda along each dimension 
        std::size_t dim = 0; 
        std::apply(
          [&](const auto&... args){
            (set_dim_boundaries_imp(args, dim++), ...);
          }, 
          m_list
        ); 
      } // end if constexpr(step)
    } // end SolBeforeStep 

    template<FDStep_Type STEP = FDStep_Type::EXPLICIT>
    void SolAfterStep(double t, const SharedConstMeshXD& mesh, StridedRef_t Sol) const 
    {
      if(mesh->numDims() != sizeof...(BCPairs_Ts)) throw std::runtime_error("BCList SolAfterStep error: MeshXD.numDims() != size of tuple / list of 1D BCs "); 
      
      if constexpr(STEP == FDStep_Type::EXPLICIT)
      {
        // this lambda takes 1 BCPair<L,R> and applies it to Sol
        // without double assigning to corners/edges of XDim space 

        auto set_dim_boundaries = [&](const auto& bc_pair, std::size_t dim){
          // Mesh1D that this bc_pair operates on  
          const auto& mesh_1dim = mesh->getMesh1D(dim); 
          // vector of eigen stride views that "look" like Discretization1Ds along mesh_1dim 
          auto views = mesh->makeOneDimViews(Sol,dim); 

          // iterate through views that look like Mesh1D
          for(std::size_t i=0; i<views.size(); i++){
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
              bc_pair.m_left.SetSolL(t, mesh_1dim, views[i]); 
              bc_pair.m_right.SetSolR(t, mesh_1dim, views[i]); 
            }
          } // end for loop through 1 dimensional views 

        }; // end set_dim_boundaries lambda 

        // apply set_dim_boundaries lambda along each dimension 
        std::size_t dim = 0; 
        std::apply(
          [&](const auto&... args){
            (set_dim_boundaries(args, dim++), ...);
          }, 
          m_list
        ); 
      } // end if constexpr(step)
    } // end SolAfterStep<>(t,m,sol) 


  private:
    // Unreachable =========================================================== 
    template<typename PAIR_T>
    MatrixStorage_t flat_stencil(double t, const SharedConstMesh1D& mesh, const PAIR_T& bc_pair) const 
    {
      // empty singuglar row of size == mesh size 
      MatrixStorage_t result(1, mesh->size()); 

      // set the row according to right boundary condition. 
      bc_pair.m_left.SetStencilL(t, mesh, result); 

      // store it into a temp vector 
      MatrixStorage_t temp(1, mesh->size()); 

      // set the row according to left boundary condition 
      bc_pair.m_right.SetStencilR(t, mesh, temp); 

      // copy left side back in to result 
      MatrixStorage_t::InnerIterator it(temp, 0); 
      for(; it; ++it) result.coeffRef(it.row(),it.col()) = it.value(); 

      return result; 
    }; 
}; // end class BCList 


} // end namespace Fds 

#endif // BCListXD.hpp 




  //   // Member Functions =======================================================

  //   // set a Matrixs' row according to list. making it an implicit stencil  
  //   virtual void SetStencil(MatrixStorage_t& Mat, const std::shared_ptr<const LinOps::MeshXD>& mesh) const override 
  //   {
  //     // check args are compaitble ---------------------------------
  //     if(list.empty()) throw std::runtime_error("list must be non empty!"); 
  //     if(list.size()!=mesh->numDims()) throw std::invalid_argument("length of boundary condition list must be == to # of numDims of MeshXD"); 
  //     // if(list.size() > 2) throw std::runtime_error("BoundaryCondXD doesn't support numDims >= 3 yet!"); 

  //     // initializations ---------------------------------------------
  //     // s1 = ith numDims size, s2 = cumulative product of all numDims < ith dim
  //     std::size_t s1=mesh->sizeOfDim(0), s2=1; 
  //     // lambda to give an NxN identity matrix 
  //     MatrixStorage_t cap_mat; 
  //     auto I = [&cap_mat](std::size_t N) -> const MatrixStorage_t& {cap_mat.resize(N,N); cap_mat.setIdentity(); return cap_mat;}; 
  //     // resulting mask for replacing rows in Mat input ---------------
  //     MatrixStorage_t mask(s1,s1);  
      
  //     // base case: 
  //     // set stencil's first/last rows maually with BcPtr_t  
  //     list[0].SetStencil(mask,mesh->getMesh1D(0)); 

  //     // recursive case: 
  //     for(std::size_t ith_dim=1; ith_dim<mesh->numDims(); ith_dim++)
  //     {
  //       // take previous mesh to higher dimension 
  //       s1 = mesh->sizeOfDim(ith_dim); 
  //       MatrixStorage_t mask_temp = Eigen::KroneckerProductSparse(I(s1), mask); 
        
  //       // fill new empty rows with ith BcPtr_t pair
  //       MatrixStorage_t sp_diag = make_SparseDiag( flat_stencil(list[ith_dim],mesh->getMesh1D(ith_dim)), s2 ); 
  //       MatrixStorage_t fill_rows = Eigen::KroneckerProductSparse(sp_diag, I(s1));
  //       fill_stencil(mask_temp, fill_rows); 

  //       // take ownership of higher dim temp with mask 
  //       mask = std::move(mask_temp); 

  //       // // increment the cumulative product 
  //       s2 *= s1; 
  //     }

  //     // overwrite rows in Mat with nonzero rows from mask 
  //     overwrite_stencil(Mat, mask); 

  //     // void return type
  //   }
  
  //   // Set time of each stored BC in 1D 
  //   virtual void SetTime(double t)
  //   {
  //     for(auto& bc : list) bc.SetTime(t); 
  //   }
  // private:
  //   // Unreachable =========================================================== 
  //   MatrixStorage_t flat_stencil(const BCPair& p, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const 
  //   {
  //     // empty singuglar row of size == mesh size 
  //     MatrixStorage_t result(1, mesh->size()); 

  //     // set the row according to left boundary condition. 
  //     p.pair.first->SetStencilL(result, mesh); 

  //     // store it into a temp vector 
  //     std::vector<double> temp(result.valuePtr(), result.valuePtr()+result.nonZeros()); 

  //     // set the row according to right boundary condition 
  //     p.pair.second->SetStencilR(result, mesh); 

  //     // copy left side back in to result 
  //     for(auto i=0; i<temp.size(); i++) result.coeffRef(0,i) = temp[i]; 

  //     return result; 
  //   }; 

/* 
// order of priority for BCListXD.list to be applied 
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

*/