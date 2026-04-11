// TimeDerivBase.hpp
//
//
//
// JAF 1/15/2026 

#ifndef TIMEDERIVBASE_H
#define TIMEDERIVBASE_H 

#include<memory>
#include "../LinOps/LinOpTraits.hpp" // fdm::Matrix
#include "TExprTraits.hpp"

namespace fdm{
  namespace texprs{
  
// Base Definition =====================================================
template<typename Derived, std::size_t N>
class TimeDerivBase
{
  public:
    // Type Defs --------------------------
    struct is_timederiv_tag{}; 
    
    // Member Data ----------------------------------
    static const std::size_t maxOrder = N; // up to nth derivative in time 

    // Constructors + Destructor ============================================
    TimeDerivBase()=default; 
    TimeDerivBase(const TimeDerivBase& other)=default; 
    ~TimeDerivBase()=default;

    // Member Funcs ========================================================== 

    // looking at container v as matrix with NCols entries per row. 
    template<std::size_t ithCol, std::size_t nCols, typename Cont>
    decltype(auto) coeffAt(const Cont& v) const 
    {
      return static_cast<const Derived*>(this)->template coeffAt<ithCol,nCols,Cont>(v); 
    } 

    // packs this pointer into iterable tuple object. note: SumExpr will override it
    auto toTuple() &
    {
      // std::cout << "Lvalue Base .toTuple()" << std::endl; 
      return std::tie(static_cast<Derived&>(*this)); // lvalue -> reference
    }
    auto toTuple() &&
    {
      // std::cout << "Rvalue Base .toTuple()" << std::endl; 
      return std::make_tuple(std::move(static_cast<Derived&&>(*this))); // rvalue -> move
    }

    // Operators ================================
    // Binary Addition Ut + Utt (Lvalue)-----------------------
    template<typename R, typename = std::enable_if_t<texprs::traits::is_timederiv_crtp<R>::value>>
    auto operator+(R&& rhs) & 
    {
      return make_SumExpr(static_cast<Derived&>(*this), std::forward<R>(rhs)); 
    }
    // (Rvalue) 
    template<typename R, typename = std::enable_if_t<texprs::traits::is_timederiv_crtp<R>::value>>
    auto operator+(R&& rhs) &&
    {
      return make_SumExpr(std::move(static_cast<Derived&>(*this)), std::forward<R>(rhs)); 
    }

    // Unary Negation Ut -> -Ut (Lvalue) -----------------------
    auto operator-() &
    {
      // delegate to Operator*() from CoeffMultExpr.hpp 
      return (-1.0) * static_cast<Derived&>(*this); 
    }
    // (Rvalue) 
    auto operator-() &&
    {
      // delegate to Operator*() from CoeffMultExpr.hpp 
      return (-1.0) * std::move(static_cast<Derived&>(*this)); 
    }

    // Binary Subtraction (Ut - Utt) -> (Ut) + (-Utt) (Lvalue) ----------------------------------
    template<typename R, typename = std::enable_if_t<texprs::traits::is_timederiv_crtp<R>::value>>
    auto operator-(R&& rhs) & 
    {
      return (*this) + (-std::forward<R>(rhs)); 
    }
    // (Rvalue)
    template<typename R, typename = std::enable_if_t<texprs::traits::is_timederiv_crtp<R>::value>>
    auto operator-(R&& rhs) && 
    {
      return std::move(*this) + (-std::forward<R>(rhs)); 
    }
}; 

  } // end namespace texprs 
} // end namespace fdm 

#endif