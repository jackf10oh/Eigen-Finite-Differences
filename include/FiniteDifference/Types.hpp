// Types.hpp 
//
// forward declarations of main fdm classes + aliases of key Eigen types 
//
// JAF 4/23/2026 

#ifndef FDM_TYPES_H
#define FDM_TYPES_H 

namespace fdm{

// forward declare ------ 
class Mesh; 

// helpful aliases ------
using Scalar = double; // might use this more consistently in the future... 
using Vector = Eigen::Matrix<fdm::Scalar, Eigen::Dynamic, 1>; 
using CSRMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>; // Compressed Sparse Row (CSR) Matrix
using DiagMatrix = Eigen::DiagonalMatrix<fdm::Scalar,Eigen::Dynamic>; // Diagonal Matrix 
using StridedRef = typename Eigen::Ref<fdm::Vector, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
using StrideView =  Eigen::Map<fdm::Vector, Eigen::Unaligned, Stride>;

} // end namespace fdm 

#endif // Types.hpp 