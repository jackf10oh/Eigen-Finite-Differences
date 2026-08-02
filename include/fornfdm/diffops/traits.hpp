// traits.hpp
//
// traits used by diffops subdirectory 
//
// JAF 4/14/2026 

#ifndef FORNFDM_DIFFOPS_TRAITS_H
#define FORNFDM_DIFFOPS_TRAITS_H

#include<cstdint>
#include<type_traits>
#include<Eigen/Core>
#include "../types.hpp"
#include "../traits.hpp"

namespace fornfdm{ 
namespace linops{
namespace internal{

// Base of all traits_impl<>. specialized by individual classes 
template<class T>
struct traits_impl{
  static constexpr bool is_linop = false; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = 0;  
  static constexpr bool is_timedep = false;
  static constexpr int direction = -1;
  static constexpr std::size_t max_order = 0;
}; 

// standard tratis of any SparseMatrix type ------------------------------- 
template<typename _Scalar, int _Options, typename _StorageIndex>
struct traits_impl< Eigen::SparseMatrix<_Scalar, _Options, _StorageIndex> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = 0; 
  static constexpr bool is_timedep = false; 
  static constexpr int direction = -1; 
  static constexpr std::size_t max_order = 0; 
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
  static constexpr std::size_t max_arity = std::max( traits_impl<std::decay_t<L>>::max_arity, traits_impl<std::decay_t<R>>::max_arity); // records maximum number of dims L/R needs to execute its callable 
  static constexpr bool is_timedep = ( traits_impl<std::decay_t<L>>::is_timedep || traits_impl<std::decay_t<R>>::is_timedep); // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = -1; // by default mixing operators direction falls back to eigen... 
  static constexpr std::size_t max_order = std::max(traits_impl<L>::max_order,traits_impl<R>::max_order); // highest order of derivative in the expression 
}; 

// Products are similar ----- 
template<class L, class R, int Options>
struct traits_impl< Eigen::Product<L,R,Options> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = std::max( traits_impl<std::decay_t<L>>::max_arity, traits_impl<std::decay_t<R>>::max_arity); // records maximum number of dims L/R needs to execute its callable 
  static constexpr bool is_timedep = ( traits_impl<std::decay_t<L>>::is_timedep || traits_impl<std::decay_t<R>>::is_timedep); // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = -1; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t max_order = std::max(traits_impl<L>::max_order,traits_impl<R>::max_order); // highest order of derivative in the expression 
}; 

// traits of Unary Expressions --------------------------------- 
template<class Op, class T>
struct traits_impl< Eigen::CwiseUnaryOp<Op, T> >
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = true; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = traits_impl<std::decay_t<T>>::max_arity; 
  static constexpr bool is_timedep = traits_impl<std::decay_t<T>>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits_impl<std::decay_t<T>>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t max_order = traits_impl<std::decay_t<T>>::max_order; // highest order of derivative in the expression 
}; 

template<class T>
using traits = traits_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

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
struct NestedStorage<T, std::void_t<decltype(fornfdm::linops::internal::traits<T>::is_linop)>>
{ 
  using traits_t = typename fornfdm::linops::internal::traits<T>; 
  typedef typename std::conditional<
    (traits_t::is_linop && !traits_t::is_unarop && !traits_t::is_binop && !traits_t::is_ternop) && (std::is_lvalue_reference<T>::value),
    T,
    typename std::remove_reference<T>::type
  >::type type; 
};

// Given two tags for NodeSelector, promote them to guarantee both minimum nodes. 
// must be specialized by actual node selectors
template<class T, class U>
struct promote_node_selector_tags
{
  constexpr static bool is_match = std::is_same_v<T,U>;  
  using type =  std::conditional_t<is_match, T, void>; 
}; 

// given std::index_sequence<> get the count
template<class T> 
struct count_orders{};

template<std::size_t... Idxs>
struct count_orders<std::index_sequence<Idxs...>>{ static constexpr std::size_t value = sizeof...(Idxs); };

// Given std::index_sequence<> get the maximum ------------------
template<std::size_t... Idxs> 
struct maximum_order_impl{};

template<std::size_t i>
struct maximum_order_impl<i>{ static constexpr std::size_t value = i; };

template<std::size_t i, std::size_t... trailing>
struct maximum_order_impl<i,trailing...>{ static constexpr std::size_t value = std::max(i, maximum_order_impl<trailing...>::value); };

template<class T>
struct maximum_order{};

template<std::size_t... Idxs>
struct maximum_order<std::index_sequence<Idxs...>>{ static constexpr std::size_t value = maximum_order_impl<Idxs...>::value; };

// Given a std::size_t, place at front of std::index_sequence<>
template<std::size_t I, class T>
struct append_order{};

template<std::size_t I, std::size_t... Idxs>
struct append_order<I, std::index_sequence<Idxs...>>{ using type = std::index_sequence<I, Idxs...>; };

// Given 2 sorted std::index_sequence<> merge them into 1
template<class LeftSeq, class RightSeq>
struct merge_orders{};

// terminating case 
template<>
struct merge_orders<std::index_sequence<>, std::index_sequence<>>{ using type = std::index_sequence<>; };

// if either side runs out of entries. append all remaining entries of other sequence
template<std::size_t... Rs>
struct merge_orders<std::index_sequence<>, std::index_sequence<Rs...>>{ using type = std::index_sequence<Rs...>; };

template<std::size_t... Ls>
struct merge_orders<std::index_sequence<Ls...>, std::index_sequence<>>{ using type = std::index_sequence<Ls...>; };

template<std::size_t L, std::size_t... Ls, std::size_t R, std::size_t... Rs>
struct merge_orders<std::index_sequence<L, Ls...>, std::index_sequence<R, Rs...>>
{
  using type = std::conditional_t<
    // if L < R, then append left entry to result sequence
    (L<R),
    typename append_order<L, typename merge_orders<std::index_sequence<Ls...>, std::index_sequence<R, Rs...>>::type>::type,
    std::conditional_t<
      // if R < L, then append left entry to result sequence
      (R<L),
      typename append_order<R, typename merge_orders<std::index_sequence<L,Ls...>, std::index_sequence<Rs...>>::type>::type,
      // if L==R, then apend exactly 1 entry to result sequence
      typename append_order<L, typename merge_orders<std::index_sequence<Ls...>, std::index_sequence<Rs...>>::type>::type
    >
  >;
};

// given an std::index_sequence<>, locate the positon of value I
template<std::size_t order, class T, class = void>
struct locate_order{};

template<std::size_t order, std::size_t I, std::size_t... Idxs>
struct locate_order<order, std::index_sequence<I, Idxs...>, std::enable_if_t<(order==I)>>
{
  static constexpr std::size_t value = 0; 
};

template<std::size_t order, std::size_t I, std::size_t... Idxs>
struct locate_order<order, std::index_sequence<I, Idxs...>, std::enable_if_t<(order!=I)>>
{
  // increment by 1, keep searching for "order"
  static constexpr std::size_t value = 1 + locate_order<order,std::index_sequence<Idxs...>>::value; 
};

// maps a derived type to its corresponding base class tag. 
// by default everything gets mapped to itself.  
template<class Derived, class = void>
struct map_to_base_tag{ using type = void; };

// if the operator uses direction zero(+ max args called <= 1), then it can be evaluated by left kronecker
struct LeftKroneckerTag{};
struct TimeDepLeftKroneckerTag : LeftKroneckerTag{}; // same tag but for time dependent case
// if the operator uses zero callable args, it can be evaluated by double kronecker
struct DoubleKroneckerTag{};
struct TimeDepDoubleKroneckerTag : DoubleKroneckerTag{}; // same tag but for time dependent case
// if the operator uses uses max arg called >= 1 (>= 2 when direction 0), we can use a weight caching evaluator
struct StoredWeightsTag{};
struct TimeDepStoredWeightsTag : StoredWeightsTag{}; // same tag but for time dependent case

} // end namespace internal

// Forward Declarations ----------
template<class Derived, typename TagType = typename internal::map_to_base_tag<Derived>::type>
class PartialDerivBase;

template<class ArgType, class = void>
class TimeEvaluation;

namespace internal{
// Given a type T detect if it is derived from PartialDerivBase<T> or is the same as PartialDerivBase<T>
template<class T>
struct is_partialderiv_crtp_helper : std::false_type{}; 

template<class T, class U>
struct is_partialderiv_crtp_helper<fornfdm::linops::PartialDerivBase<T,U>> : std::true_type{}; 

template<class T>
struct is_partialderiv_crtp_impl : std::disjunction<std::is_base_of<fornfdm::linops::PartialDerivBase<T>,T>, is_partialderiv_crtp_helper<T>>{}; 

template<class T>
using is_partialderiv_crtp = is_partialderiv_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 


} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm

#endif // traits.hpp (diffops)