// functors.hpp 
//
// Function Object declarations for NwiseUnaryOp and NwiseBinaryOp 
//
// JAF 4/23/2026 

#ifndef FORNFDM_DIFFOPS_FUNCTORS_H
#define FORNFDM_DIFFOPS_FUNCTORS_H 

#include<type_traits>
#include<Eigen/Core> // need the declarations from BinaryFunctors.h in Core 

namespace fornfdm{
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

// utility to convert any Eigen CwiseBinaryOp / CwiseUnaryOp's functor to whatever created the expression
template<typename T>
struct ConvertedFO: public T
{
  ConvertedFO(const T& functor)
    : T(functor)
  {}
};

// scalar_opposite_op ---
template<typename ScalarType>
struct ConvertedFO<Eigen::internal::scalar_opposite_op<ScalarType>> : public UnaryNegateFO
{
  ConvertedFO(const Eigen::internal::scalar_opposite_op<ScalarType>& functor)
    : UnaryNegateFO()
  {}
};

// scalar_sum_op ---
template<typename LhsScalar,typename RhsScalar>
struct ConvertedFO<Eigen::internal::scalar_sum_op<LhsScalar, RhsScalar>> : public BinaryAdditionFO
{
  ConvertedFO(const Eigen::internal::scalar_sum_op<LhsScalar, RhsScalar>& functor)
    : BinaryAdditionFO()
  {}
};

// scalar_difference_op ---
template<typename LhsScalar,typename RhsScalar>
struct ConvertedFO<Eigen::internal::scalar_difference_op<LhsScalar, RhsScalar>> : public BinarySubtractionFO
{
  ConvertedFO(const Eigen::internal::scalar_difference_op<LhsScalar, RhsScalar>& functor)
    : BinarySubtractionFO()
  {}
};

// scalar_product_op ---
template<typename LhsScalar,typename RhsScalar>
struct ConvertedFO<Eigen::internal::scalar_product_op<LhsScalar, RhsScalar>> : public BinaryMultiplyFO
{
  ConvertedFO(const Eigen::internal::scalar_product_op<LhsScalar, RhsScalar>& functor)
    : BinaryMultiplyFO()
  {}
};

// bind2nd_op --- 
template<typename BinaryOp>
struct ConvertedFO<Eigen::internal::bind2nd_op<BinaryOp>> : public ConvertedFO<BinaryOp>
{
  typename BinaryOp::second_argument_type m_arg; 
  ConvertedFO(const Eigen::internal::bind2nd_op<BinaryOp>& functor)
    : ConvertedFO<BinaryOp>(functor), m_arg(functor.m_value)
  {}
  template<class LhsType>
  auto operator()(const LhsType& lhs) const
  {
    return ConvertedFO<BinaryOp>::operator()(lhs, m_arg);
  }
};

// bind1st_op --- 
template<typename BinaryOp>
struct ConvertedFO<Eigen::internal::bind1st_op<BinaryOp>> : public ConvertedFO<BinaryOp>
{
  typename BinaryOp::first_argument_type m_arg; 
  ConvertedFO(const Eigen::internal::bind1st_op<BinaryOp>& functor)
    : ConvertedFO<BinaryOp>(functor), m_arg(functor.m_value)
  {}
  template<class RhsType>
  auto operator()(const RhsType& rhs) const
  {
    return ConvertedFO<BinaryOp>::operator()(m_arg, rhs);
  }
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#endif // functors.hpp 