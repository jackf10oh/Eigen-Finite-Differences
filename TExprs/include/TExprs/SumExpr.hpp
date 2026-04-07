// SumExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSSUMEXPR_H
#define LHSSUMEXPR_H 

#include "TExprTraits.hpp"
#include "TimeDerivBase.hpp"

namespace TExprs{

// ====================================================== 
template<typename... Args>
class SumExpr : public TimeDerivBase<SumExpr<Args...>, TExprs::internal::variadicFoldMaximum<Args::maxOrder...>::value> 
{
  public:
    std::tuple<Args...> m_args; 
  public:
    // Constructors + Destructor ================================
    SumExpr()=delete; 
    SumExpr(std::tuple<Args...> tup_init, std::size_t order) 
      : m_args( tup_init ) 
    {} 
    // destructor 
    ~SumExpr()=default; 

    // Member Funcs ====================================================
    template<std::size_t ithCol, std::size_t nCols, typename Cont>
    decltype(auto) coeffAt(const Cont& v) const = delete; 

    // Lvalue toTuple 
    auto& toTuple() & 
    { return m_args; }

    // Rvalue toTuple -> move args 
    auto toTuple() &&
    { return std::move(m_args); }

    std::string toString() const {return "hi from sum"; }; 
};

template<typename LHS, typename RHS>
auto make_SumExpr(LHS&& lhs, RHS&& rhs)
{
  std::size_t m = std::max(lhs.Order(), rhs.Order()); 
  auto left_tup = std::forward<LHS>(lhs).toTuple(); 
  auto right_tup = std::forward<RHS>(rhs).toTuple(); 
  auto cat = std::tuple_cat(left_tup, right_tup); 
  return SumExpr(cat, m); 
}

} // end namespace TExprs 

#endif // SumExpr.hpp 