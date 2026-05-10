// FillStencil.hpp
//
//
//
// JAF 1/2/2025 

#ifndef FDM_UTILS_FILLSTENCIL_H
#define FDM_UTILS_FILLSTENCIL_H

#include<Eigen/Sparse> 

namespace fdm{
  namespace utils{

/* template<typename Scalar, typename Index, typename Derived> 
void fill_stencil(Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& A, const Eigen::SparseMatrixBase<Derived>& mask)
{
  // check A and mask are same size 
  static_assert(Derived::IsRowMajor, "Error in fill_stencil: Mask expression must be row major!"); 
  if(A.cols()!=mask.cols() || A.rows()!=mask.rows()) throw std::invalid_argument("Error in fill_stencil: rows/cols of stencil / mask must be =="); 
  
  // loop throw rows (cols) of mask 
  for(auto i=0; i < mask.rows(); ++i)
  {
    bool A_empty = true; 
    typename Eigen::InnerIterator it(A,i);
    // if any entries are non zero in A's row (col)  
    for(; it; ++it){
      if(it.value()!=0.0){
        A_empty=false; 
        break; 
      }
    }; 

    // if A has no entries in row i 
    if(A_empty){ 
      Eigen::InnerIterator A_it(A,i);
      for(; A_it; ++A_it){
        A_it.valueRef() = Scalar(0.0); 
      }

      // fill A's row (col) with entries in mask's row (col) 
      typename Eigen::InnerIterator M_it(mask.derived(), i); 
      for(; M_it; ++M_it){
        A.coeffRef(M_it.row(),M_it.col()) = M_it.value(); 
      }
    }
    // else: just leave row i in A untouched 
  }
};

*/ 

template<typename Scalar, typename Index, typename Derived> 
void overwrite_stencil(Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& A, const Eigen::SparseMatrixBase<Derived>& mask)
{
  // check A and mask are same size 
  static_assert(Derived::IsRowMajor, "Error in fill_stencil: Mask expression must be row major!"); 
  if(A.cols()!=mask.cols() || A.rows()!=mask.rows()) throw std::invalid_argument("Error in fill_stencil: rows/cols of stencil / mask must be =="); 

  // loop throw rows (cols) of mask 
  for(auto ith_row=0; ith_row < A.rows(); ++ith_row)
  {    
    // !!! Assume that mask has no Explicit Zeros !!! 
    Eigen::InnerIterator<Derived> M_it(mask.derived(),ith_row); 
    // if mask has any non zero entries  
    if(M_it) { 
      // set all values of A's row to zero 
      A.uncompress(); 
      A.innerNonZeroPtr()[ith_row] = 0;

      // fill A's row (col) with entries in mask's row (col) 
      for(; M_it; ++M_it){
        A.coeffRef(M_it.row(),M_it.col()) = M_it.value(); 
      }
    }
  }
};

  } // end namespace utils 
} // end namespace fdm 

#endif // FillStencil.hpp