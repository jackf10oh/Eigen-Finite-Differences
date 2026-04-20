// Traits.hpp
//
//
//
// JAF 4/14/2026 

#ifndef FDMTRAITS_H
#define FDMTRAITS_H

#include<cstdint>
#include<string>
#include<type_traits> // decay_t 
#include<complex>
#include<Eigen/Core>
#include<Eigen/src/Core/util/Macros.h>
#include<Eigen/src/Core/util/Constants.H> 
#include<Eigen/src/Core/util/ForwardDeclarations.h>  // CwiseUnaryOp, CwiseBinaryOp, 
#include<Eigen/src/Core/EigenBase.h> 
#include<Eigen/src/SparseCore/SparseUtil.h> // forward declares SparseMatrix<...> 
#include<Eigen/src/SparseCore/CompressedStorage.h>
#include<Eigen/src/SparseCore/SparseCompressedBase.h>
// class CompressedStorage; 

// #include<Eigen/SparseCore> can't include before plugin macro takes effect! 

namespace fdm{ 
// forward declare ------ 
class Mesh; 

// helpful aliases ------
using Scalar = double; // might use this more consistently in the future... 
using RowMajorMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;
using CSRMatrix = RowMajorMatrix; // Column Sparse Row (CSR) Matrix
using StridedRef = typename Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
using StrideView =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride>;

namespace internal{

// template trait to detect if T<> and U<> are the same templates 
template<template<typename...> class T, template<typename...> class U>
struct is_same_template : std::false_type{}; 

template<template<typename...> class T>
struct is_same_template<T,T> : std::true_type{}; 

template<class T>
struct traits_impl{}; 

// standard tratis of any SparseMatrix type ------------------------------- 
template<typename _Scalar, int _Options, typename _StorageIndex>
struct traits_impl< Eigen::SparseMatrix<_Scalar, _Options, _StorageIndex> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = 0; 
  static constexpr bool is_timedep = false; 
  static constexpr int direction = -1; 
  static constexpr std::size_t maxOrder = 0; 
}; 

// standard traits for derived from Eigen::SparseMatrixBase or SparseCompressedBase ------------------------------------
template<class Derived>
struct traits_impl< Eigen::SparseMatrixBase<Derived> > : traits_impl< std::decay_t<Derived> >{}; 

template<class Derived>
struct traits_impl< Eigen::SparseCompressedBase<Derived> > : traits_impl< std::decay_t<Derived> >{}; 

// traits of a  Binary Expression ------------------------------------
template<class Op, class L, class R>
struct traits_impl< Eigen::CwiseBinaryOp<Op,L,R> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = std::max( traits_impl<std::decay_t<L>>::max_num_args_called, traits_impl<std::decay_t<R>>::max_num_args_called); // records maximum number of dims L/R needs to execute its callable 
  static constexpr bool is_timedep = ( traits_impl<std::decay_t<L>>::is_timedep || traits_impl<std::decay_t<R>>::is_timedep); // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = -1; // by default mixing operators direction falls back to eigen... 
  static constexpr std::size_t maxOrder = std::max(traits_impl<L>::maxOrder,traits_impl<R>::maxOrder); // highest order of derivative in the expression 
}; 

// Products are similar ----- 
template<class L, class R, int Options>
struct traits_impl< Eigen::Product<L,R,Options> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = std::max( traits_impl<std::decay_t<L>>::max_num_args_called, traits_impl<std::decay_t<R>>::max_num_args_called); // records maximum number of dims L/R needs to execute its callable 
  static constexpr bool is_timedep = ( traits_impl<std::decay_t<L>>::is_timedep || traits_impl<std::decay_t<R>>::is_timedep); // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = -1; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = std::max(traits_impl<L>::maxOrder,traits_impl<R>::maxOrder); // highest order of derivative in the expression 
}; 

// traits of Unary Expressions --------------------------------- 
template<class Op, class T>
struct traits_impl< Eigen::CwiseUnaryOp<Op, T> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = true; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = traits_impl<std::decay_t<T>>::max_num_args_called; 
  static constexpr bool is_timedep = traits_impl<std::decay_t<T>>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits_impl<std::decay_t<T>>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = traits_impl<std::decay_t<T>>::maxOrder; // highest order of derivative in the expression 
}; 

template<class T>
using traits = traits_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

} // end namespace internal 
} // end namespace fdm

#endif // Traits.hpp 