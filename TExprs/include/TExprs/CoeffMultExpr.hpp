// CoeffMultExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSCOEFFMULTEXPR_H
#define LHSCOEFFMULTEXPR_H

#include "TimeDerivBase.hpp"
#include "NthTimeDeriv.hpp" 
#include<LinOps/LinOpTraits.hpp> // LinOps::traits::Storage<> 

namespace TExprs{
namespace internal{

// ======================================================
template<typename COEFF_T, typename RHS_T>
class CoeffMultExpr : public TimeDerivBase<CoeffMultExpr<COEFF_T,RHS_T>> 
{
  public:
    // Type Defs --------------------- 
    using Lhs_t = typename LinOps::traits::Storage<COEFF_T>::type; 
    using Rhs_t = typename std::conditional_t<
      std::conjunction_v<
        std::is_same<std::remove_cv_t<std::remove_reference_t<RHS_T>>, NthTimeDeriv>, 
        std::is_lvalue_reference<RHS_T>
      >, 
      RHS_T, 
      std::remove_reference_t<RHS_T>
    >; 

    // Member Data ---------------
    Lhs_t m_coeff; 
    Rhs_t m_rhs; 
    /* store the right hand side somehow... 
    if its an lvalue NthTimeDeriv store a reference. 
    if its a rvalue NthTimeDeriv just store a copy.
    All other types (CoeffMultExpr) are stored by Copy*/

  public:
    // Constructors + Destructor ====================================
    // CoeffMultExpr()=delete; // unnecessary?  
    CoeffMultExpr(Lhs_t c_init, Rhs_t rhs_init)
      : m_coeff(c_init), m_rhs(rhs_init), TimeDerivBase<CoeffMultExpr>(rhs_init.Order())
    {}
    CoeffMultExpr(const CoeffMultExpr& other)=default; 
    // destructor 
    ~CoeffMultExpr()=default; 

    // Member Funcs =================================================== 
    // L/R Getters
    Lhs_t& Lhs(){return m_coeff; } 
    Rhs_t& Rhs(){return m_rhs; } 
    const Lhs_t& Lhs() const {return m_coeff; } 
    const Rhs_t& Rhs() const {return m_rhs; } 

    // Must be implemented for TimeDerivBase.  Get coeff * rhs's CoeffAt() value  
    template<typename Cont>
    auto CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      if constexpr(LinOps::traits::is_coeffop_crtp<Lhs_t>::value){
        return m_coeff.getMat() * m_rhs.CoeffAt(v,n_nodes_per_row,ith_node);  
      }
      else if constexpr(std::is_same<double, std::remove_cv_t<std::remove_reference_t<Lhs_t>>>::value){
        return m_coeff * m_rhs.CoeffAt(v,n_nodes_per_row,ith_node); 
      }
      else{
        throw std::runtime_error("LhsCoeffMultExpr error: COEFF_T is not a double or CoeffOp."); 
      }
    } 

    // set_mesh overrides TimeDerivBase 
    template<typename ANYMESHPTR_T>
    void set_mesh(ANYMESHPTR_T m)
    {
      // is m_coeff is a LinOp and not a double call its set_mesh(); 
      if constexpr(LinOps::traits::is_coeffop_crtp<Lhs_t>::value){
        m_coeff.set_mesh(m); 
      }
      // m_rhs is either a CoeffMultExpr of NthTimeDeriv... 
      m_rhs.set_mesh(m); 
    } // end set_mesh(m) 

    // set_time(t) overrides TimeDerivBase
    void setTime(double t)
    {
      if constexpr(LinOps::traits::is_coeffop_crtp<Lhs_t>::value)
      {
        m_coeff.setTime(t); 
      }
      m_rhs.setTime(t); 
    }

    // using LhsBase<LhsCoeffMultExpr<COEFF_T>>::toTuple; 
    std::string toString() const {return "hi from coeffMult!"; }; 
}; // end class CoeffMultExpr

} // end namespace internal 
} // end namespace TExprs 

template<
  typename Lhs, 
  typename Rhs
>
std::enable_if_t<
  std::conjunction_v<
    LinOps::traits::is_coeffop_crtp<Lhs>,
    TExprs::traits::is_timederiv_crtp<Rhs>
  >, 
  TExprs::internal::CoeffMultExpr<Lhs,Rhs>
> operator*(Lhs&& c, Rhs&& rhs)
{
  // false if Rhs is any form of SumExpr. We don't want to mess with expressions like c*(A+B)
  static_assert(std::tuple_size<decltype(rhs.toTuple())>::value == 1, "operator*(c,TimeDeriv) only meant for single TimeDeriv"); 
  return TExprs::internal::CoeffMultExpr<Lhs,Rhs>(std::forward<Lhs>(c), std::forward<Rhs>(rhs)); 
} 

// Operator for scalar c, TimeDeriv Ut making expression c*Ut 
template<
  typename SCALAR_T, 
  typename TDERIV_RHS
>
std::enable_if_t<
  std::conjunction_v<
    TExprs::traits::is_timederiv_crtp<TDERIV_RHS>,
    std::is_arithmetic<std::remove_reference_t<std::remove_cv_t<SCALAR_T>>>
  >, 
  TExprs::internal::CoeffMultExpr<SCALAR_T,TDERIV_RHS>
> operator*(SCALAR_T&& c, TDERIV_RHS&& rhs)
{
  return TExprs::internal::CoeffMultExpr<SCALAR_T,TDERIV_RHS>(std::forward<SCALAR_T>(c), std::forward<TDERIV_RHS>(rhs)); 
}

#endif // LhsCoeffMultExpr.hpp 