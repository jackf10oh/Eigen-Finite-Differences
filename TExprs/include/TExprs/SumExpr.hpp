// SumExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSSUMEXPR_H
#define LHSSUMEXPR_H 

#include "TimeDerivBase.hpp"

namespace TExprs{
namespace internal{

// ====================================================== 
template<typename... Args>
class SumExpr : public TimeDerivBase<SumExpr<Args...>> 
{
  public:
    std::tuple<Args...> m_args; 
  public:
    // Constructors + Destructor ================================
    SumExpr()=delete; 
    SumExpr(std::tuple<Args...> tup_init, std::size_t order) 
      : m_args( tup_init ), TimeDerivBase<SumExpr<Args...>>(order) 
    {} 
    // destructor 
    ~SumExpr()=default; 

    // Member Funcs ====================================================
    template<typename Cont>
    auto CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const = delete;
    template<typename ANYMESHPTR_T>
    void set_mesh(ANYMESHPTR_T m)=delete; 
    void setTime(double t)=delete; 
    auto& toTuple() &
    {
      // std::cout << "Lvalue SumExpr .toTuple()" << std::endl; 
      return m_args; // lvalue -> reference
    }
    auto toTuple() &&
    {
      // std::cout << "Rvalue SumExpr .toTuple()" << std::endl; 
      return std::move(m_args); // rvalue -> rvalue ref
    }
    std::string toString() const {return "hi from sum"; }; 
};

template<typename LHS, typename RHS>
auto make_sumexpr_helper(LHS&& lhs, RHS&& rhs)
{
  std::size_t m = std::max(lhs.Order(), rhs.Order()); 
  auto left_tup = std::forward<LHS>(lhs).toTuple(); 
  auto right_tup = std::forward<RHS>(rhs).toTuple(); 
  auto cat = std::tuple_cat(left_tup, right_tup); 
  return SumExpr(cat, m); 
}

} // end namespace internal
} // end namespace TExprs 

#endif // SumExpr.hpp 