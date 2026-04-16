// EigenPlugin.hpp 
// 
// header file that adds functionality to 
// ALL eigen sparse matrices to be able to 
//
// 1. ) take in a shared_ptr to a mesh discretization of a spatial domain 
// 
// 2. ) take in a double representing the current time 
// 
// JAF 4/13/2026 

// #ifndef FDM_EIGEN_SPARSE_BASE_PLUGIN // no include guards. shouldn't be directly included anywhere else... 
// #define FDM_EIGEN_SPARSE_BASE_PLUGIN

// using fdm::traits; 
// class fdm::Mesh; 

public:
// Type Defs ------------- 
struct is_linop_tag{}; 
// struct is_time_dep{}; // put this in derived classes if they override setTime() in a meaniful way. i.e. update state of linear operator

// Member Functions ================================================================== 
void setMesh(const std::shared_ptr<const fdm::Mesh>& m) 
{
  using D = std::remove_reference_t<std::remove_cv_t<Derived>>; 
  // TODO ternary operators? 
  if constexpr(fdm::traits::is_binop<Derived>::value){
    // binary expressions hook lhs/rhs
    std::cout <<"binop SetMesh! "; 
    if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value){
      std::cout <<"L hooked, "; 
      derived().lhs().const_cast_derived().setMesh(m);      
    }
    if constexpr(fdm::traits::is_linop<typename Derived::Rhs>::value){
      std::cout <<"R Hooked" << std::endl;
      std::cout << "RHS type: " << typeid(decltype(derived().rhs())).name() << std::endl; 
      std::cout << "RHS  const cast type: " << typeid(decltype(derived().rhs().const_cast_derived())).name() << std::endl; 
      derived().rhs().const_cast_derived().setMesh(m);
    }
    std::cout << std::endl;
  }
  else if constexpr(fdm::traits::is_unarop<Derived>::value){
    // unary expressions hook nestedExpression 
    std::cout <<"unarop SetMesh!"; 
    if constexpr(fdm::traits::is_linop<typename Derived::XprTypeNested>::value){
      std::cout << "nested hooked"; 
      derived().nestedExpression().const_cast_derived().setMesh(m); 
    }
    std::cout << std::endl; 
  }
  else{
    // leaf matrices resize... 
    // std::size_t s = m->sizesProduct(); 
    // const_cast<D&>(derived()).resize(s,s); 
  }
}

std::shared_ptr<const fdm::Mesh> getMesh() const 
{
  // TODO ternary operators? 
  if constexpr(fdm::traits::is_binop<Derived>::value){
    // binary expressions hook lhs/rhs
    if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value && fdm::traits::is_linop<typename Derived::Rhs>::value){
      auto result = derived().lhs().derived().getMesh(); 
      return result ? result : derived().rhs().derived().getMesh();  
    }
    else if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value){
      return derived().lhs().derived().getMesh();
    }
    else if constexpr(fdm::traits::is_linop<typename Derived::Rhs>::value){
      return derived().rhs().derived().getMesh();
    }       
    else{
      return nullptr; 
    }
  }
  else if constexpr(fdm::traits::is_unarop<Derived>::value){
    // unary expressions hook nestedExpression 
    if constexpr(fdm::traits::is_linop<typename Derived::XprTypeNested>::value){
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

void setTime(double t)
{
  using D = std::remove_reference_t<std::remove_cv_t<Derived>>; 
  // TODO ternary operators? 
  if constexpr(fdm::traits::is_binop<Derived>::value){
    // binary expressions hook lhs/rhs
    if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value){
      derived().lhs().const_cast_derived().setTime(t);      
    }
    if constexpr(fdm::traits::is_linop<typename Derived::Rhs>::value){
      derived().rhs().const_cast_derived().setTime(t);      
    }       
  }
  else if constexpr(fdm::traits::is_unarop<Derived>::value){
    // unary expressions hook nestedExpression 
    if constexpr(fdm::traits::is_linop<typename Derived::XprTypeNested>::value){
      derived().nestedExpression().const_cast_derived().setTime(t);      
    }
  }
  // else leaf matrices do nothing by default; 
}

double getTime() const 
{
  // TODO ternary operators? 
  if constexpr(fdm::traits::is_binop<Derived>::value){
    // binary expressions hook lhs/rhs
    if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value && fdm::traits::is_linop<typename Derived::Rhs>::value){
      auto result = derived().lhs().getTime(); 
      return (result != -1.0) ? result : derived().rhs().getTime();  
    }
    else if constexpr(fdm::traits::is_linop<typename Derived::Lhs>::value){
      return derived().lhs()().getTime();
    }
    else if constexpr(fdm::traits::is_linop<typename Derived::Rhs>::value){
      return derived().rhs().getTime();
    }       
    else{
      return -1.0; 
    }
  }
  else if constexpr(fdm::traits::is_unarop<Derived>::value){
    // unary expressions hook nestedExpression 
    if constexpr(fdm::traits::is_linop<typename Derived::XprTypeNested>::value){
      return derived().nestedExpression().getTime(); 
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

// Doesn't work...
// static constexpr bool isTimeDep = is_time_dep<decltype(std::declval<Derived>().derived())>::value; 

// goes into SparseMatrix<> class

// #endif // EigenPlugin.hpp 