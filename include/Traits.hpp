// Traits.hpp
//
//
//
// JAF 4/14/2026 

#ifndef FDMTRAITS_H
#define FDMTRAITS_H

#include<cstdint>
#include<string>
#include<type_traits>
#include<complex>
#include<Eigen/Core>
#include<Eigen/src/Core/util/Macros.h>
#include<Eigen/src/Core/util/Constants.H> 
#include<Eigen/src/Core/util/ForwardDeclarations.h>  // CwiseUnaryOp, CwiseBinaryOp, Eigen::StorageOptions::RowMajor
#include<Eigen/src/SparseCore/SparseUtil.h> // forward declares SparseMatrix<...> 
// #include<Eigen/src/SparseCore/CompressedStorage.h>
// #include<Eigen/src/SparseCore/SparseCompressedBase.h>
// class CompressedStorage; 

// #include<Eigen/SparseCore> can't include before plugin macro takes effect! 

namespace fdm{ 
// forward declare ------ 
class Mesh; 

// helpful aliases ------
using Matrix = Eigen::SparseMatrix<double,Eigen::StorageOptions::RowMajor>; 
using StridedRef = typename Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
using StrideView =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride>;

namespace traits{
// Unary Ops --------------------------------------------------------------
template<typename T>
struct is_unarop_impl : std::false_type{}; 

template<typename Op, typename T>
struct is_unarop_impl<Eigen::CwiseUnaryOp<Op,T>> : std::true_type{};

template<typename T>
using is_unarop = is_unarop_impl<T>; 

// Binary Ops --------------------------------------------------------------
template<typename T>
struct is_binop_impl : std::false_type{}; 

template<typename Op, typename L, typename R>
struct is_binop_impl<Eigen::CwiseBinaryOp<Op, L, R>> : std::true_type{}; 

template<typename L, typename R, int Option>
struct is_binop_impl<Eigen::Product<L, R, Option>> : std::true_type{}; 

template<typename T>
using is_binop = is_binop_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

// Ternary Ops --------------------------------------------------------------
template<typename T>
struct is_ternarop_impl : std::false_type{}; 

template<typename Op, typename A, typename B, typename C>
struct is_ternarop_impl<Eigen::CwiseTernaryOp<Op,A,B,C>> : std::true_type{};

template<typename T>
using is_ternarop = is_ternarop_impl<T>;

// detect if T is an Eigen sparse matrix --------------------------------- 
template<typename T, typename = void>
struct is_linop_impl : std::false_type{}; 

template<typename T>
struct is_linop_impl<T, std::void_t<typename T::is_linop_tag> > : std::true_type{}; 

// TODO Diagonal matrices are linops 

template<typename T>
using is_linop = is_linop_impl< std::remove_reference_t<std::remove_cv_t<T>> >; 


// Time Dependencies -------------------------------------------------------------- 
// TODO ternary ops?
template<typename T, typename = void>
struct is_time_dep_impl : std::false_type{}; 

template<typename T>
struct is_time_dep_impl<T, std::void_t<typename T::is_timedep_tag> > : std::true_type{}; 

template<typename Op, typename T>
struct is_time_dep_impl<Eigen::CwiseUnaryOp<Op,T>> : is_time_dep_impl<typename Eigen::CwiseUnaryOp<Op,T>::XprTypeNested>{};

template<typename Op, typename L, typename R>
struct is_time_dep_impl<Eigen::CwiseBinaryOp<Op,L,R>> : std::disjunction<
  is_time_dep_impl<typename Eigen::CwiseBinaryOp<Op,L,R>::Lhs>,
  is_time_dep_impl<typename Eigen::CwiseBinaryOp<Op,L,R>::Rhs>
>{};

template<typename L, typename R, int Option>
struct is_time_dep_impl<Eigen::Product<L,R,Option>> : std::disjunction<
  is_time_dep_impl<typename Eigen::Product<L,R,Option>::Lhs>,
  is_time_dep_impl<typename Eigen::Product<L,R,Option>::Rhs>
>{};

template<typename T, typename = void>
struct is_time_dep : std::false_type{}; 

template<typename T>
struct is_time_dep<T, std::enable_if_t<is_linop<T>::value> > : is_time_dep_impl<std::remove_reference_t<std::remove_cv_t<T>>>{}; 

} // end namespace traits 
} // end namespace fdm

#endif // Traits.hpp 