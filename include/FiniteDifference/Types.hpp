// Types.hpp 
//
// forward declarations of main fornfdm classes + aliases of key Eigen types 
//
// JAF 4/23/2026 

#ifndef FORNFDM_TYPES_H
#define FORNFDM_TYPES_H 

#include<cstdint>
#include<string>
#include<type_traits> // decay_t 
#include<complex>
#include<Eigen/Core>
#include<Eigen/src/Core/NumTraits.h> // convert scalar -> real 
#include<Eigen/src/Core/util/Macros.h>
#include<Eigen/src/Core/util/Constants.H> 
#include<Eigen/src/Core/util/ForwardDeclarations.h>  // CwiseUnaryOp, CwiseBinaryOp, 
#include<Eigen/src/Core/EigenBase.h> 
#include<Eigen/src/SparseCore/SparseUtil.h> // forward declares SparseMatrix<...> 
#include<Eigen/src/SparseCore/CompressedStorage.h>
#include<Eigen/src/SparseCore/SparseCompressedBase.h>

namespace fornfdm{

// forward declare ------ 
class Mesh; 
template<std::size_t N> class Coordinate; 

// helpful aliases ------
#ifndef FORNFDM_CUSTOM_SCALAR
using Scalar = double; // might use this more consistently in the future...
#else
using Scalar = FORNFDM_CUSTOM_SCALAR;
#endif 

#ifndef FORNFDM_CUSTOM_REAL
using Real = typename Eigen::NumTraits<fornfdm::Scalar>::Real;  
#else 
using Real = FORNFDM_CUSTOM_REAL; 
#endif
 
using Vector = typename Eigen::Matrix<fornfdm::Scalar, Eigen::Dynamic, 1>; 
using CSRMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>; // Compressed Sparse Row (CSR) Matrix
using DiagMatrix = Eigen::DiagonalMatrix<fornfdm::Scalar,Eigen::Dynamic>; // Diagonal Matrix 
using StridedRef = typename Eigen::Ref<fornfdm::Vector, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
using StrideView =  Eigen::Map<fornfdm::Vector, Eigen::Unaligned, Stride>;

} // end namespace fornfdm 

#endif // Types.hpp 