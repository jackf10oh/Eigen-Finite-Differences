// plugin.hpp 
// 
// header file that adds functionality to 
// ALL eigen sparse matrices to be able to 
//
// 1a. ) take in a shared_ptr to a mesh discretization of a spatial domain 
// 1b. ) return a shared_ptr<const Mesh> to the last mesh set 
// 
// 2a. ) take in a double representing the current time 
// 2b. ) return the last double set as time  
//
// 3 ) cast to const Derived&. 
// -- Used by fornfdm expression to explicitly use Eigen operators. i.e. +,-,* over fornfdm operators 
// 
// JAF 4/13/2026 

#ifndef EIGEN_SPARSEMATRIXBASE_H

#ifndef FORNFDM_PLUGIN_H
#define FORNFDM_PLUGIN_H

// fornfdm plugin dependencies
#include "diffops/traits.hpp"

// set Eigen's plugin as this file  
#ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <fornfdm/plugin.hpp>
#endif

#endif // FORNFDM_PLUGIN_H

#else // EIGEN_SPARSEMATRIXBASE_H

#ifndef FORNFDM_DIFFOPS_TRAITS_H
#error "<fornfdm/plugin.hpp> depends on <fornfdm/diffops/traits.hpp>"
#endif

public: 
// Member Functions ================================================================== 
const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }

void setMesh(const std::shared_ptr<const fornfdm::Mesh>& m) 
{
  if constexpr(fornfdm::linops::internal::traits<Derived>::is_binop){
    // binary expressions hook lhs/rhs
    if constexpr(fornfdm::linops::internal::traits<typename Derived::Lhs>::is_linop){
      derived().lhs().const_cast_derived().setMesh(m);      
    }
    if constexpr(fornfdm::linops::internal::traits<typename Derived::Rhs>::is_linop){
      derived().rhs().const_cast_derived().setMesh(m);
    }
  }
  else if constexpr(fornfdm::linops::internal::traits<Derived>::is_unarop){
    // unary expressions hook nestedExpression 
    if constexpr(fornfdm::linops::internal::traits<typename Derived::XprTypeNested>::is_linop){
      derived().nestedExpression().const_cast_derived().setMesh(m); 
    }
  }
  else{
    // leaf matrices resize?  
    // std::size_t s = m->sizesProduct(); 
    // const_cast<D&>(derived()).resize(s,s); 
  }
}

std::shared_ptr<const fornfdm::Mesh> getMesh() const 
{
  if constexpr(fornfdm::linops::internal::traits<Derived>::is_binop){
    // binary expressions hook lhs/rhs
    constexpr bool left_returns = fornfdm::linops::internal::traits<typename Derived::Lhs>::is_linop; 
    constexpr bool right_returns = fornfdm::linops::internal::traits<typename Derived::Rhs>::is_linop; 
    if constexpr(left_returns && right_returns){
      auto result = derived().lhs().derived().getMesh(); 
      return (result != nullptr) ? result : derived().rhs().derived().getMesh();  
    }
    else if constexpr(left_returns){
      return derived().lhs().derived().getMesh();
    }
    else if constexpr(right_returns){
      return derived().rhs().derived().getMesh();
    }       
    else{
      return nullptr; 
    }
  }
  else if constexpr(fornfdm::linops::internal::traits<Derived>::is_unarop){
    // unary expressions hook nestedExpression 
    if constexpr(fornfdm::linops::internal::traits<typename Derived::XprTypeNested>::is_linop){
      return derived().nestedExpression().derived().getMesh(); 
    }
    else{
      return nullptr; 
    }
  }
  else{
    // leaf matrices return nullptr by default 
    return nullptr;  
  }
}

void setTime(fornfdm::Real t)
{
  if constexpr(fornfdm::linops::internal::traits<Derived>::is_binop){
    // binary expressions hook lhs/rhs
    if constexpr(fornfdm::linops::internal::traits<typename Derived::Lhs>::is_linop){
      derived().lhs().const_cast_derived().setTime(t);      
    }
    if constexpr(fornfdm::linops::internal::traits<typename Derived::Rhs>::is_linop){
      derived().rhs().const_cast_derived().setTime(t);      
    }       
  }
  else if constexpr(fornfdm::linops::internal::traits<Derived>::is_unarop){
    // unary expressions hook nestedExpression 
    if constexpr(fornfdm::linops::internal::traits<typename Derived::XprTypeNested>::is_unarop){
      derived().nestedExpression().const_cast_derived().setTime(t);      
    }
  }
  // else leaf matrices do nothing by default; 
}

fornfdm::Real getTime() const 
{
  if constexpr(fornfdm::linops::internal::traits<Derived>::is_binop){
    // binary expressions hook lhs/rhs
    constexpr bool left_returns = fornfdm::linops::internal::traits<typename Derived::Lhs>::is_linop; 
    constexpr bool right_returns = fornfdm::linops::internal::traits<typename Derived::Rhs>::is_linop; 
    if constexpr(left_returns && right_returns){
      auto result = derived().lhs().derived().getTime(); 
      return (result != -1.0) ? result : derived().rhs().derived().getTime();  
    }
    else if constexpr(left_returns){
      return derived().lhs().derived().getTime();
    }
    else if constexpr(right_returns){
      return derived().rhs().derived().getTime();
    }       
    else{
      return -1.0; 
    }
  }
  else if constexpr(fornfdm::linops::internal::traits<Derived>::is_unarop){
    // unary expressions hook nestedExpression 
    if constexpr(fornfdm::linops::internal::traits<typename Derived::XprTypeNested>::is_linop){
      return derived().nestedExpression().derived().getTime(); 
    }
    else{
      return -1.0; 
    }
  }
  else{
    // leaf matrices return -1.0 by default 
    return -1.0;  
  }
}

#endif // EIGEN_SPARSEMATRIXBASE_H

// plugin.hpp 