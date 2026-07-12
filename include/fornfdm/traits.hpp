// traits.hpp
//
// traits used by all subdirectory of fornfdm library 
//
// JAF 4/14/2026 

#ifndef FORNFDM_TRAITS_H
#define FORNFDM_TRAITS_H

#include<cstdint>
#include<type_traits> // decay_t 
#include "types.hpp"

namespace fornfdm{ 
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
    std::is_convertible<ReturnType, fornfdm::Scalar>, 
    std::is_convertible<fornfdm::Scalar, Args>
    ...
  >; 
}; 

// operator() const 
template<class ClassType, class ReturnType, typename... Args>
struct callable_traits<ReturnType(ClassType::*)(Args...) const>
{
  static constexpr std::size_t arity = sizeof...(Args); 
  static constexpr bool maps_scalars_to_scalar = std::conjunction_v<
    std::is_convertible<ReturnType, fornfdm::Scalar>, 
    std::is_convertible<fornfdm::Scalar, Args>
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

} // end namespace internal 
} // end namespace fornfdm

#endif // traits.hpp 