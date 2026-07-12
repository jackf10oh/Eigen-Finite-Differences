// traits.hpp
//
// compile time traits for TExpr types 
//
// JAF 4/4/2026 

#ifndef FORNFDM_TEXPR_TRAITS_H
#define FORNFDM_TEXPR_TRAITS_H

#include<cstdint>
#include<tuple>
#include<type_traits>
#include "../types.hpp"

namespace fornfdm{
  namespace texprs{

template<typename Derived, std::size_t N>
class TimeDerivBase; 

template<std::size_t N>
class NthTimeDeriv; 

template<typename Coeff, typename TimeDeriv>
class CoeffMultExpr; 

// traits =======================================================================================
// See if a given type T derives from TimeDerivBase ---------------------------------------------
namespace internal{
template<typename T, typename = void> 
struct is_timederiv_crtp_impl : public std::false_type{}; 

template<typename T> 
struct is_timederiv_crtp_impl<T, std::void_t<typename T::is_timederiv_tag>>: public std::true_type{}; 
} // end namespace internal 

namespace traits{
template<typename T>
using is_timederiv_crtp = texprs::internal::is_timederiv_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 
} // end namespace traits

// See if a given type T derives from CoeffMultExpr ------------------------------
namespace internal{
template<typename T> 
struct is_coeffmult_crtp_impl : public std::false_type{}; 

template<typename L, typename R> 
struct is_coeffmult_crtp_impl<CoeffMultExpr<L,R>>: public std::true_type{}; 
} // end namespace internal 

namespace traits{
template<typename T>
using is_coeffmult_crtp = texprs::internal::is_coeffmult_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 
} // end namespace traits

// Given a type T, determine how to store it in CoeffMultExpr ------------------------------------ 
namespace traits{ 

template<typename T, typename = void>
struct Storage
{ 
  using type = T; 
}; 

template<typename T>
struct Storage<T, std::enable_if_t<is_timederiv_crtp<T>::value>>
{
  using type = std::conditional_t<
    std::is_lvalue_reference<T>::value && !(is_coeffmult_crtp<T>::value), 
    T,
    std::remove_reference_t<T>
  >; 
}; 

} // end namespace traits 

// Given size_t i,j,... get the maximum -------------------------------------------- 
namespace internal{
  
template<std::size_t... is>
struct variadicFoldMaximum{/* no value for empty args*/}; 

template<std::size_t i>
struct variadicFoldMaximum<i>
{
  static constexpr std::size_t value = i; 
}; 

template<std::size_t i, std::size_t... trailing>
struct variadicFoldMaximum<i,trailing...>
{
  static constexpr std::size_t value = std::max(i, variadicFoldMaximum<trailing...>::value); 
};

} // end namespace internal 

// given a type T. detect if it returns a fornfdm::Real from coeffAt<>() or some other matrix expression. 
namespace traits{

template<typename TIMEDERIV_T>
struct TimeDerivTraits
{
  using _coeffatreturntype_unclean_ = decltype(std::declval<TIMEDERIV_T>().template coeffAt<0,0>(std::declval<std::array<fornfdm::Real,0>>()));  
  using CoeffAtReturnType = std::remove_cv_t<std::remove_reference_t<_coeffatreturntype_unclean_>>;
}; 

template<typename T>
struct coeffat_returns_real_impl : std::false_type{};

template<std::size_t nthOrder>
struct coeffat_returns_real_impl<texprs::NthTimeDeriv<nthOrder>> : std::true_type{}; 

template<class Lhs, class Rhs>
struct coeffat_returns_real_impl<texprs::CoeffMultExpr<Lhs,Rhs>> : std::conjunction<
  std::is_convertible<decltype(std::declval<Lhs>() * fornfdm::Real()),fornfdm::Real>, 
  coeffat_returns_real_impl<std::remove_cv_t<std::remove_reference_t<Rhs>>>
>{}; 

template<class T>
using coeffat_returns_real = coeffat_returns_real_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 

template<typename T>
struct coeffat_returns_other : public std::negation<coeffat_returns_real<T>>{}; 

} // end namespace traits

// Given a tuple <A,B,C,...> return a subset where PRED<X>::value returns true.   
namespace traits{ 

template<template<typename...> class PRED, typename TUP_T>
auto filter_tup(TUP_T&& tup)
{
  auto filter_lam = [](auto&& elem){
    using T = decltype(elem); 
    using CLEAN_T = std::remove_reference_t<std::remove_cv_t<T>>;
    if constexpr(PRED<CLEAN_T>::value){
      // singleton tuple ( elem )
      if constexpr(std::is_lvalue_reference<T>::value){
        return std::tie(elem); 
      }
      else{
        return std::make_tuple(std::move(elem)); 
      }
    }
    else{
      // empty tuple ( _ )
      return std::tuple<>{}; 
    }
  }; 

  auto singletons = std::apply(
    [&](auto&&... elems){
      return std::make_tuple(filter_lam(std::forward<decltype(elems)>(elems))...); 
    }, 
    tup
  ); 

  return std::apply(
    [](auto&&... elems){
      return std::tuple_cat(std::forward<decltype(elems)>(elems)...); 
    }, 
    singletons
  ); 
}

} // end namespace traits 

  } // end namespace texprs 
} // end namespace fornfdm 

#endif 