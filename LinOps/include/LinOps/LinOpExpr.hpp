// LinOpExpr.hpp
//
// Complie time expression template that stores Lhs, Rhs, BinOp (UnarOp), and Mesh
//
// JAF 12/7/2025

#ifndef LINOPEXPR_H
#define LINOPEXPR_H

#include<type_traits>
#include "LinOpBase.hpp"
#include "LinOpTraits.hpp"

namespace linops{

// Expression of L,R, BinOp, + Mesh =====================================================
template<typename L, typename R, typename BinaryOp>
class LinOpExpr : public LinOpMixIn< LinOpExpr<L, R, BinaryOp> >, public LinOpBase1D< LinOpExpr<L, R, BinaryOp> >, public LinOpBaseXD< LinOpExpr<L, R, BinaryOp> >
{
  public:
    // Type Defs --------------------------------------
    using LStorage = typename traits::Storage<L>::type;
    using RStorage = typename traits::Storage<R>::type;
    using Operator = BinaryOp; 

    // shadows LinOpMixin's isTimeDep. 
    static constexpr bool isTimeDep = traits::is_time_dep<L>::value || traits::is_time_dep<R>::value;  

  private:
    // Member Data ---------------------------------------------
    LStorage m_lhs;
    RStorage m_rhs;
    const BinaryOp m_binop; 
    
  public:
    // Constructors + Destructor =============================================================
    
    // default 
    LinOpExpr()=delete; 

    // from Lhs, Rhs, BinOp, + Mesh 
    LinOpExpr(LStorage A, RStorage B)
      : m_lhs(A), m_rhs(B), m_binop()
    {
      constexpr bool both_linop = traits::is_linop_crtp<L>::value && traits::is_linop_crtp<R>::value;
      constexpr bool both_1d = traits::is_1dim_linop_crtp<L>::value && traits::is_1dim_linop_crtp<R>::value;
      constexpr bool both_xd = traits::is_xdim_linop_crtp<L>::value && traits::is_xdim_linop_crtp<R>::value;
      if constexpr(both_linop){
        static_assert(both_1d || both_xd, "Error constructing LinOpExpr: tried to mix strictly 1D operator with strictly XD operator"); 
      }
    };

    // copy 
    LinOpExpr(const LinOpExpr& other)=default; 

    // destructor
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================

    // returns combination bin_op(A,B) of 2 stored linops ----------------
    decltype(auto) asMatrix() const
    {
      return m_binop(m_lhs, m_rhs); 
    };

    // fixes ambiguous .apply() .set_mesh()? 
    using LinOpBase1D<LinOpExpr<L,R,BinaryOp> >::apply; 
    using LinOpBaseXD<LinOpExpr<L,R,BinaryOp> >::apply; 
    using LinOpBase1D<LinOpExpr<L,R,BinaryOp>>::setMesh; 
    using LinOpBaseXD<LinOpExpr<L,R,BinaryOp>>::setMesh; 

    // Expr Only ==============================================================

    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage& getLhs(){ return m_lhs; }; 
    const LStorage& getLhs() const { return m_lhs; }; 
    // rhs 
    RStorage& getRhs(){ return m_rhs; }; 
    const RStorage& getRhs() const { return m_rhs; }; 
   
};

// Specialization for Unary operators --------------------------------------------
template<typename L, typename UnaryOp>
class LinOpExpr<L, void, UnaryOp> : public LinOpMixIn< LinOpExpr<L, void, UnaryOp> >, public LinOpBase1D< LinOpExpr<L, void, UnaryOp> >, public LinOpBaseXD< LinOpExpr<L, void, UnaryOp> >
{
  public:
    // Type Defs ------------------------------------------------------------------
    using LStorage = typename traits::Storage<L>::type;
    using RStorage = void; // not storing a second argument anymore 
    using Operator = UnaryOp; 

  private:
    // Member Data -------------------------------------------------------------
    LStorage m_lhs;
    // RStorage m_rhs; // not storing a second argument anymore 
    const UnaryOp m_unarop; 
    
  public:
    // Constructors / Destructors ===============================================
    // default 
    LinOpExpr()=delete; 
    // from lhs, unar_op, + mesh
    LinOpExpr(LStorage A)
      : m_lhs(A), m_unarop()
    {
      constexpr bool A_is_linop = traits::is_linop_crtp<L>::value;
      static_assert(A_is_linop, "Error constructing Unary LinOpExpr: single linop A is not linop!"); 
    };
    // copy 
    LinOpExpr(const LinOpExpr& other)=default; 
    // destructors
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================

    // returns Op( A ) of 1 stored linops --------------------------------------
    decltype(auto) asMatrix() const
    {
      // return m_binop(m_lhs, m_rhs); 
      return m_unarop(m_lhs); 
    };

    // fixes ambiguous .apply() ?
    using LinOpBase1D<LinOpExpr<L,void,UnaryOp>>::apply; 
    using LinOpBaseXD<LinOpExpr<L,void,UnaryOp>>::apply;  

    // Expr Only ==============================================================

    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage& getLhs(){ return m_lhs; };     
    const LStorage& getLhs() const { return m_lhs; };     
    // since nothing is stored. but still declared so that decltype(getRhs()) is still usable
    void getRhs() const {};    
    
};

} // end namespace linops 

#endif // LinOpExpr.hpp