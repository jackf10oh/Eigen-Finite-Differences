// CoeffMultExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSCOEFFMULTEXPR_H
#define LHSCOEFFMULTEXPR_H

#include<LinOps/LinOpTraits.hpp> // LinOps::traits::Storage_t<>, is_linop_crtp<> 
#include "TimeDerivBase.hpp"
#include "TExprTraits.hpp"

namespace TExprs{

// ======================================================
template<typename Coeff, typename TimeDeriv>
class CoeffMultExpr : public TimeDerivBase<CoeffMultExpr<Coeff,TimeDeriv>, std::remove_reference_t<std::remove_cv_t<TimeDeriv>>::maxOrder> 
{
  public:
    // Type Defs --------------------- 
    /* if its an lvalue TimeDeriv store a reference. 
    if its a rvalue NthTimeDeriv just store a copy. */ 
    using LStorage = typename LinOps::traits::Storage_t<Coeff>::type; 
    using RStorage = typename TExprs::traits::Storage<TimeDeriv>::type; 

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
      if constexpr(LinOps::traits::is_coeffop_crtp<LStorage>::value){
        return m_coeff.GetMat() * m_rhs.template coeffAt<ithCol,nCols,Cont>(v);  
      }
      else if constexpr(std::is_same<double, std::remove_cv_t<std::remove_reference_t<LStorage>>>::value){
        return m_coeff * m_rhs.template coeffAt<ithCol,nCols,Cont>(v); 
      }
      else{
        throw std::runtime_error("CoeffMultExpr error: LHS (Coeff) is not a double or CoeffOp."); 
      }
    } 

    /* // set_mesh overrides TimeDerivBase 
    template<typename AnyMesh>
    void set_mesh(const std::shared_ptr<AnyMesh>& m)
    {
      // is m_coeff is a LinOp and not a double call its set_mesh(); 
      if constexpr(LinOps::traits::is_coeffop_crtp<LStorage>::value){
        m_coeff.set_mesh(m); 
      }
      // m_rhs is either a CoeffMultExpr of NthTimeDeriv... 
      m_rhs.set_mesh(m); 
    } // end set_mesh(m) 

    // set_time(t) overrides TimeDerivBase
    void SetTime(double t)
    {
      if constexpr(LinOps::traits::is_coeffop_crtp<LStorage>::value)
      {
        m_coeff.SetTime(t); 
      }
      m_rhs.SetTime(t); 
    }

    */ 

    // using LhsBase<LhsCoeffMultExpr<COEFF_T>>::toTuple; 
    std::string toString() const {return "hi from coeffMult!"; }; 
}; // end class CoeffMultExpr

// Operator for CoeffOp c, TimeDeriv Ut making expression c*Ut 
template<
  typename Lhs, 
  typename Rhs
>
std::enable_if_t<
  std::conjunction_v<
    LinOps::traits::is_coeffop_crtp<Lhs>,
    TExprs::traits::is_timederiv_crtp<Rhs>
  >, 
  TExprs::CoeffMultExpr<Lhs,Rhs>
> operator*(Lhs&& c, Rhs&& rhs)
{
  // false if Rhs is any form of SumExpr. We don't want to mess with expressions like c*(A+B)
  static_assert(std::tuple_size<decltype(rhs.toTuple())>::value == 1, "operator*(c,TimeDeriv) only meant for single TimeDeriv"); 
  return TExprs::CoeffMultExpr<Lhs,Rhs>(std::forward<Lhs>(c), std::forward<Rhs>(rhs)); 
} 

// Operator for Scalar c, TimeDeriv Ut making expression c*Ut 
template<
  typename Scalar, 
  typename TimeDeriv, 
  typename = std::enable_if_t<
    std::conjunction_v<
      TExprs::traits::is_timederiv_crtp<TimeDeriv>,
      std::is_arithmetic<std::remove_reference_t<std::remove_cv_t<Scalar>>>
    >
  > 
> 
auto operator*(Scalar&& c, TimeDeriv&& rhs)
{
  return TExprs::CoeffMultExpr<Scalar,TimeDeriv>(std::forward<Scalar>(c), std::forward<TimeDeriv>(rhs)); 
}

} // end namespace TExprs 

#endif // LhsCoeffMultExpr.hpp 