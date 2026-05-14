// Traits.hpp
//
// traits used by all subdirectory of fdm library 
//
// JAF 4/14/2026 

#ifndef FDM_TRAITS_H
#define FDM_TRAITS_H

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
// #include<Eigen/SparseCore> can't include before plugin macro takes effect! 
#include "Types.hpp"

namespace fdm{ 
namespace internal{

// traits around a callable type F ---------------------------------------- 
template<typename T> 
struct callable_traits : public callable_traits<decltype(&T::operator())>{}; 

// operator() -- non const
template<class ClassType, class ReturnType, typename... Args>
struct callable_traits<ReturnType(ClassType::*)(Args...)>
{
  static constexpr std::size_t arity = sizeof...(Args); 
  static constexpr bool maps_scalars_to_scalar = std::conjunction_v<
    std::is_convertible<ReturnType, fdm::Scalar>, 
    std::is_convertible<fdm::Scalar, Args>
    ...
  >; 
}; 

// operator() const 
template<class ClassType, class ReturnType, typename... Args>
struct callable_traits<ReturnType(ClassType::*)(Args...) const>
{
  static constexpr std::size_t arity = sizeof...(Args); 
  static constexpr bool maps_scalars_to_scalar = std::conjunction_v<
    std::is_convertible<ReturnType, fdm::Scalar>, 
    std::is_convertible<fdm::Scalar, Args>
    ...
  >; 
}; 

template<typename T>
struct BindFirst : BindFirst<decltype(&T::operator())>
{
  // Constructor
  using BindFirst<decltype(&T::operator())>::BindFirst; 
}; 

template<class ClassType, class ReturnType, typename First, typename... Trailing>
class BindFirst<ReturnType(ClassType::*)(First, Trailing...) const>
{
  private:
    // Type Defs ----------------------- 
    using MemFn = ReturnType(ClassType::*)(First, Trailing...) const;  
    
    // Member Data ------------------------------ 
    const MemFn mem_fn = &ClassType::operator(); 
    const std::remove_reference_t<ClassType> binded; 
  public:
    std::remove_reference_t<First> captured_arg;

    // Constructor ==================================
    BindFirst(ClassType c, First x, MemFn f = &ClassType::operator())
      : binded(std::move(c)), captured_arg(std::move(x)), mem_fn(f)
    {} 

    // Member Functions / Operators ---------------- 
    ReturnType operator()(Trailing... args) const 
    {
      return (binded.*mem_fn)(captured_arg, args...); 
    }
}; 

// detect if two classes are the same. irregardless of non class template params 
template<class A, class B>
struct is_matching_impl : std::is_same<A,B>{}; 

// template<class A>
// struct is_matching_impl<A,A> : std::true_type{}; 

template<template< auto... > class T, auto... Ls, auto... Rs>
struct is_matching_impl<T<Ls...>, T<Rs...>> : std::true_type{};

template<class A, class B>
using is_matching = is_matching_impl<std::remove_cv_t<std::remove_reference_t<A>>, std::remove_cv_t<std::remove_reference_t<B>>>; 

} // end namespace internal 
} // end namespace fdm

#endif // Traits.hpp 