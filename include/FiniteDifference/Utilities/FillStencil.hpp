// FillStencil.hpp
//
//
//
// JAF 1/2/2025 

#ifndef FILLSTENCIL_H
#define FILLSTENCIL_H 

#include<Eigen/Sparse> 

template<typename Scalar, typename Index> 
void fill_stencil(Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& A, const Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& mask)
{
  // check A and mask are same size 
  if(A.outerSize()!=mask.outerSize() || A.innerSize()!=mask.innerSize()) throw std::invalid_argument("Outer Size and Inner Size of Matrix A and Mask Matrix must be ==."); 
  
  // loop throw rows (cols) of mask 
  for(auto i=0; i < mask.outerSize(); i++)
  {
    typedef typename Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>::InnerIterator InnerIterator;
    
    bool A_empty = true; 

    // if any entries are non zero in A's row (col)  
    for(InnerIterator it(A,i); it; it.operator++()) if(it.value()!=0.0) A_empty=false; 

    // if A has any non zero entries do nothing to the row (col) 
    if(A_empty) { 
      // else: set all values of A's row (col) to zero 
      for(InnerIterator A_it(A,i); A_it; A_it.operator++()) A_it.valueRef() = Scalar(0.0); 

      // fill A's row (col) with entries in mask's row (col) 
      for(InnerIterator it(mask,i); it; it.operator++()){
        A.coeffRef(it.row(),it.col()) = it.value(); 
      }
    }
  }
};

template<typename Scalar, typename Index> 
void overwrite_stencil(Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& A, const Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>& mask)
{
  // check A and mask are same size 
  if(A.outerSize()!=mask.outerSize() || A.innerSize()!=mask.innerSize()) throw std::invalid_argument("Outer Size and Inner Size of Matrix A and Mask Matrix must be ==."); 
  
  // loop throw rows (cols) of mask 
  for(auto i=0; i < mask.outerSize(); i++)
  {
    typedef typename Eigen::SparseMatrix<Scalar,Eigen::RowMajor,Index>::InnerIterator InnerIterator;
    
    bool mask_empty = true; 

    // if any entries are non zero in A's row (col)  
    for(InnerIterator it(mask,i); it; it.operator++()) if(it.value()!=0.0) mask_empty=false; 

    // if mask has any non zero entries  
    if(!mask_empty) { 
      // set all values of A's row (col) to zero 
      for(InnerIterator A_it(A,i); A_it; A_it.operator++()) A_it.valueRef() = Scalar(0.0); 

      // fill A's row (col) with entries in mask's row (col) 
      for(InnerIterator it(mask,i); it; it.operator++()){
        A.coeffRef(it.row(),it.col()) = it.value(); 
      }
    }
  }
};

#endif // FillStencil.hpp