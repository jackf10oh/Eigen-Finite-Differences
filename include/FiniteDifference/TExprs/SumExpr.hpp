// SumExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSSUMEXPR_H
#define LHSSUMEXPR_H 

#include "TExprTraits.hpp"
#include "TimeDerivBase.hpp"

namespace fdm{
  namespace texprs{

// ====================================================== 
template<typename... Args>
class SumExpr : public TimeDerivBase<SumExpr<Args...>, texprs::internal::variadicFoldMaximum<std::remove_reference_t<std::remove_cv_t<Args>>::maxOrder...>::value> 
{
  public:
    std::tuple<Args...> m_args; 
  public:
    // Constructors + Destructor ================================
    SumExpr()=delete; 

    SumExpr(std::tuple<Args...> tup_init) 
      : m_args( tup_init ) 
    {} 
    
    SumExpr(const SumExpr& other)=default;

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

};

template<typename LHS, typename RHS>
auto make_SumExpr(LHS&& lhs, RHS&& rhs)
{
  static_assert(
    texprs::traits::is_timederiv_crtp<LHS>::value && texprs::traits::is_timederiv_crtp<RHS>::value,
    "SumExpr must be constructed from 2 Time Derivative Expressions!"
  );
  return SumExpr(std::tuple_cat(std::forward<LHS>(lhs).toTuple(),std::forward<RHS>(rhs).toTuple())); 
}

  } // end namespace texprs 
} // end namespace fdm 

#endif // SumExpr.hpp 