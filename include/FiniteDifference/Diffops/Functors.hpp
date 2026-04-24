// Functors.hpp 
//
// Function Object declarations for NwiseUnaryOp and NwiseBinaryOp 
//
// JAF 4/23/2026 

#ifndef FDM_DIFFOP_FUNCTORS_H
#define FDM_DIFFOP_FUNCTORS_H 

namespace fdm{
namespace linops{
namespace internal{

struct UnaryNegateFO
{
  template<typename X>
  auto operator()(const X& x) const { return -x; } 
}; 

template<typename C>
struct UnaryScalarMultiplyFO
{
  using ScalarTypeNested =
    std::conditional_t<
      std::is_lvalue_reference_v<C>,
      C,
      std::remove_reference_t<C>
    >;

  ScalarTypeNested m_scalar;

  // Constructor
  UnaryScalarMultiplyFO(C&& c)
    : m_scalar(std::forward<C>(c))
  {}

  template<typename X>
  auto operator()(const X& x) const
  {
    return m_scalar * x;
  }
};

struct BinaryAdditionFO
{
  template<typename X, typename Y>
  auto operator()(const X& x, const Y& y) const { return x + y; }
}; 

struct BinarySubtractionFO
{
  template<typename X, typename Y>
  auto operator()(const X& x, const Y& y) const { return x - y; }
}; 

struct BinaryMultiplyFO
{
  template<typename X, typename Y>
  auto operator()(const X& x, const Y& y) const { return x * y; }
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm 

#endif // Functors.hpp 