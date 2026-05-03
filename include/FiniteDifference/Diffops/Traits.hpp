// Traits.hpp
//
// Traits used by linops sub directory 
//
// JAF 4/14/2026 

#ifndef FDM_LINOP_TRAITS_H
#define FDM_LINOP_TRAITS_H

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

#include "../Types.hpp"
#include "../Traits.hpp"

// #include<Eigen/SparseCore> can't include before plugin macro takes effect! 

namespace fdm{ 
namespace linops{

// Forward Declarations ------------------------ 
template<class Derived>
struct PartialDerivBase; 

namespace internal{

// Base of all traits_impl<>. specialized by individual classes 
template<class T>
struct traits_impl{
  static constexpr bool is_linop = false; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = 0;  
  static constexpr bool is_timedep = false;
  static constexpr int direction = -1;
  static constexpr std::size_t maxOrder = 0;
}; 

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

// Given a type T detect if it is derived from PartialDerivBase<T> or is the same as PartialDerivBase<T>
template<class T>
struct is_partialderiv_crtp_helper : std::false_type{}; 

template<class T>
struct is_partialderiv_crtp_helper<fdm::linops::PartialDerivBase<T>> : std::true_type{}; 

template<class T>
struct is_partialderiv_crtp_impl : std::disjunction<std::is_base_of<fdm::linops::PartialDerivBase<T>,T>, is_partialderiv_crtp_helper<T>>{}; 

template<class T>
using is_partialderiv_crtp = is_partialderiv_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 

// Determine what type to use to nest any partial deriv 
template<class T, typename = void>
struct NestedStorage
{
  typedef typename std::conditional<
    (std::is_lvalue_reference<T>::value),
    T,
    typename std::remove_reference<T>::type
  >::type type; 
};

template<class T>
struct NestedStorage<T, std::void_t<decltype(Eigen::internal::traits<std::remove_cv_t<std::remove_reference_t<T>>>::Flags)>>
{
  typedef typename std::conditional<
    (Eigen::internal::traits<std::remove_cv_t<std::remove_reference_t<T>>>::Flags & Eigen::NestByRefBit) && (std::is_lvalue_reference<T>::value),
    T,
    typename std::remove_reference<T>::type
  >::type type; 
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm

#endif // Traits.hpp 