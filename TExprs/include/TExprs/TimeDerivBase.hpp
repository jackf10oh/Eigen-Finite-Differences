// TimeDerivBase.hpp
//
//
//
// JAF 1/15/2026 

#ifndef TIMEDERIVBASE_H
#define TIMEDERIVBASE_H 

#include<LinOps/LinearOpBase.hpp> // MatrixStorage_t

namespace TExprs{

using MatrixStorage_t = LinOps::MatrixStorage_t; 

// traits =====================================================
namespace internal{
template<typename T, typename = void> 
struct is_timederiv_crtp_impl : public std::false_type{}; 

template<typename T> 
struct is_timederiv_crtp_impl<T, std::void_t<typename T::is_timederiv_tag>>: public std::true_type{}; 
} // end namespace internal 

namespace traits{
template<typename T>
using is_timederiv_crtp = TExprs::internal::is_timederiv_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 
} // end namespace traits

namespace internal{
  
// Base Definition =====================================================
template<typename Derived>
class TimeDerivBase
{
  public:
    // Type Defs --------------------------
    struct is_timederiv_tag{}; 
  public:
    // Member Data ----------------------------------
    std::size_t m_order; // nth derivative in time 
  public:
    // Constructors + Destructor ============================================
    TimeDerivBase()=default; 
    TimeDerivBase(std::size_t order_init)
      : m_order(order_init)
    {}
    TimeDerivBase(const TimeDerivBase& other)=default; 
    ~TimeDerivBase()=default;
    // Member Funcs ========================================================== 
    // Current Order of expression  
    std::size_t Order() const {return m_order; }; 

    template<typename Cont>
    decltype(auto) CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return static_cast<const Derived*>(this)->CoeffAt(v, n_nodes_per_row, ith_node); 
    } 

    // going to need set time and set mesh so that coeff ops can interac with Lhs... 
    // default of set_mesh(m) just do nothing. TimeDeriv doesn't depend on spatial
    template<typename ANYMESHPTR_T>
    void set_mesh(ANYMESHPTR_T m){}; 

    // default of setTime() just do nothing. 
    void setTime(double t){}; 

    // packs this pointer into iterable tuple object. note: SumExpr will override it
    auto toTuple() &
    {
      // std::cout << "Lvalue Base .toTuple()" << std::endl; 
      return std::tie(static_cast<Derived&>(*this)); // lvalue -> reference
    }
    auto toTuple() &&
    {
      // std::cout << "Rvalue Base .toTuple()" << std::endl; 
      return std::make_tuple(static_cast<Derived&&>(*this)); // rvalue -> rvalue ref
    }

    std::string toString() const {return "hi from base"; }; 

    // Operators ================================
    // Binary Addition Ut + Utt (Lvalue)-----------------------
    template<typename RHS, typename = std::enable_if_t<TExprs::traits::is_timederiv_crtp<RHS>::value>>
    auto operator+(RHS&& rhs) & 
    {
      return make_sumexpr_helper(static_cast<Derived&>(*this), std::forward<RHS>(rhs)); 
    }
    // (Rvalue) 
    template<typename RHS, typename = std::enable_if_t<TExprs::traits::is_timederiv_crtp<RHS>::value>>
    auto operator+(RHS&& rhs) &&
    {
      return make_sumexpr_helper(std::move(static_cast<Derived&>(*this)), std::forward<RHS>(rhs)); 
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
    template<typename RHS, typename = std::enable_if_t<TExprs::traits::is_timederiv_crtp<RHS>::value>>
    auto operator-(RHS&& rhs)
    {
      return this->operator+(-std::forward<RHS>(rhs)); 
    }
}; 

} // end namespace internal 

} // end namespace TExprs

#endif