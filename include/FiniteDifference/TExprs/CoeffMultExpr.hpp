// CoeffMultExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSCOEFFMULTEXPR_H
#define LHSCOEFFMULTEXPR_H

#include "../LinOps/LinOpTraits.hpp" // linops::traits::Storage_t<>, is_linop_crtp<>
#include "TimeDerivBase.hpp"
#include "TExprTraits.hpp"

namespace fdm{
  namespace texprs{

// ======================================================
template<typename Coeff, typename TimeDeriv>
class CoeffMultExpr : public TimeDerivBase<CoeffMultExpr<Coeff,TimeDeriv>, std::remove_reference_t<std::remove_cv_t<TimeDeriv>>::maxOrder> 
{
  public:
    // Type Defs --------------------- 
    /* if its an lvalue TimeDeriv store a reference. 
    if its a rvalue NthTimeDeriv just store a copy. */ 
    using LStorage = typename linops::traits::Storage<Coeff>::type; 
    using RStorage = typename texprs::traits::Storage<TimeDeriv>::type; 

    // Member Data ---------------
    LStorage m_coeff; 
    RStorage m_rhs; 

  public:
    // Constructors + Destructor ====================================
    CoeffMultExpr()=delete;  
    
    CoeffMultExpr(LStorage c_init, RStorage rhs_init)
      : m_coeff(c_init), m_rhs(rhs_init)
    {}

    CoeffMultExpr(const CoeffMultExpr& other)=default; 

    // destructor 
    ~CoeffMultExpr()=default; 

    // Member Funcs =================================================== 
    // L/R Getters
    LStorage& getLhs(){return m_coeff; } 
    RStorage& getRhs(){return m_rhs; } 
    const LStorage& getLhs() const {return m_coeff; } 
    const RStorage& getRhs() const {return m_rhs; } 

    // Must be implemented for TimeDerivBase.  Get coeff * rhs's CoeffAt() value  
    template<std::size_t ithCol, std::size_t nCols, typename Cont>
    decltype(auto) coeffAt(const Cont& v) const 
    {
      if constexpr(linops::traits::is_coeffop_crtp<LStorage>::value){
        return m_coeff.asMatrix() * m_rhs.template coeffAt<ithCol,nCols,Cont>(v);  
      }
      else if constexpr(std::is_same<double, std::remove_cv_t<std::remove_reference_t<LStorage>>>::value){
        return m_coeff * m_rhs.template coeffAt<ithCol,nCols,Cont>(v); 
      }
      else{
        throw std::runtime_error("CoeffMultExpr error: LHS (Coeff) is not a double or CoeffOp."); 
      }
    } 
}; // end class CoeffMultExpr

// Operator for CoeffOp c, TimeDeriv Ut making expression c*Ut 
template<
  typename Lhs, 
  typename Rhs
>
std::enable_if_t<
  std::conjunction_v<
    linops::traits::is_coeffop_crtp<Lhs>,
    texprs::traits::is_timederiv_crtp<Rhs>
  >, 
  texprs::CoeffMultExpr<Lhs,Rhs>
> operator*(Lhs&& c, Rhs&& rhs)
{
  // false if Rhs is any form of SumExpr. We don't want to mess with expressions like c*(A+B)
  static_assert(std::tuple_size<decltype(rhs.toTuple())>::value == 1, "operator*(c,TimeDeriv) only meant for single TimeDeriv"); 
  return CoeffMultExpr<Lhs,Rhs>(std::forward<Lhs>(c), std::forward<Rhs>(rhs)); 
} 

// Operator for Scalar c, TimeDeriv Ut making expression c*Ut 
template<
  typename Scalar, 
  typename TimeDeriv, 
  typename = std::enable_if_t<
    std::conjunction_v<
      texprs::traits::is_timederiv_crtp<TimeDeriv>,
      std::is_arithmetic<std::remove_reference_t<std::remove_cv_t<Scalar>>>
    >
  > 
> 
auto operator*(Scalar&& c, TimeDeriv&& rhs)
{
  return CoeffMultExpr<Scalar,TimeDeriv>(std::forward<Scalar>(c), std::forward<TimeDeriv>(rhs)); 
}

  } // end namespace texprs 
} // end namespace fdm 

#endif // CoeffMultExpr.hpp 