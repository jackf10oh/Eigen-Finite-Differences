// CoeffMultExpr.hpp
//
// binary expressions like c*Ut 
//
// JAF 1/15/2026 

#ifndef FDM_TEXPRS_COEFFMULTEXPR_H
#define FDM_TEXPRS_COEFFMULTEXPR_H

#include "TimeDerivBase.hpp"
#include "TExprTraits.hpp"
#include "../Diffops/Traits.hpp"

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
    using LStorage = typename fdm::linops::internal::NestedStorage<Coeff>::type; 
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
      return (m_rhs.template coeffAt<ithCol,nCols,Cont>(v)) * m_coeff;  
      // return m_rhs.template coeffAt<ithCol,nCols,Cont>(v);  
    } 

}; // end class CoeffMultExpr

// Operator for coeff c, TimeDeriv Ut making expression c*Ut 
template<typename Lhs, typename Rhs, typename = std::enable_if_t<texprs::traits::is_timederiv_crtp<Rhs>::value>>
auto operator*(Lhs&& c, Rhs&& rhs)
{
  // false if Rhs is any form of SumExpr. We don't want to mess with expressions like c*(A+B)
  static_assert(std::tuple_size<decltype(rhs.toTuple())>::value == 1, "operator*(c,TimeDeriv) only meant for single TimeDeriv"); 
  std::cout << "Product Made! lhs is lval? " << std::is_lvalue_reference<Lhs>::value << " rhs is lval? " << std::is_lvalue_reference<Rhs>::value << std::endl; 
  return CoeffMultExpr<Lhs,Rhs>(std::forward<Lhs>(c), std::forward<Rhs>(rhs)); 
} 

} // end namespace texprs 
} // end namespace fdm 

#endif // CoeffMultExpr.hpp 