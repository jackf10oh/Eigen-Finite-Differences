// NwiseUnaryOp.hpp
//
//
//
// JAF 4/22/2026 

#ifndef FORNFDM_DIFFOPS_NWISEUNARYOP_H
#define FORNFDM_DIFFOPS_NWISEUNARYOP_H

#include<cstdint>
#include "../types.hpp"
#include "traits.hpp"
#include "EvaluatorBase.hpp"
#include "PartialDerivBase.hpp"
#include "EigenEvaluatorImpl.hpp"
#include "functors.hpp"

namespace fornfdm{
namespace linops{

// Forward Declaration
template<class UnaryOp, class XprType>
struct NwiseUnaryOp; 

namespace internal{

// Evaluator 
template<class UnaryOp, class XprType>
struct Evaluator<fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>> : public EvaluatorBase<fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>>
{
  using UnarXprType = fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>;  
  const UnarXprType& m_xpr; 
  Evaluator<typename UnarXprType::NestedExpression> m_nested_eval; 
  Evaluator(const fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>& xpr) : m_xpr(xpr), m_nested_eval(xpr.nestedExpression()){} 

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return [f = m_xpr.functor(), nested = m_nested_eval.createReader(coord,t)](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride)
    {
      return f(nested(weights,idx,stride));
    };
  }

  template<std::size_t N>
  class ExactReader
  {
    using NestedReader = decltype(std::declval<const Evaluator<typename UnarXprType::NestedExpression>&>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>(), std::declval<fornfdm::Real>()));
    const UnarXprType& m_xpr;
    NestedReader m_nested; 
    ExactReader(const UnarXprType& xpr, NestedReader nested) 
    : m_xpr(xpr), m_nested(nested)
    {};
    template<std::size_t... orders>
    operator()(const fornfdm::Scalar* data, std::size_t idx, std::size_t stride, std::index_sequence<orders...>) const
    {
      return m_xpr.functor()(m_nested(data, idx, stride, std::index_sequence<orders...>{}));
    }
  };

  template<std::size_t N>
  auto createExactReader(const fornfdm::Coordinate<N> coord, fornfdm::Real t) const
  {
    return ExactReader<N>(m_xpr, m_nested_eval.template createExactReader<N>(coord, t));
  }
};

// traits
template<class UnaryOp, class XprType>
struct traits_impl<fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>> : public traits_impl<XprType>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = true; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = traits<XprType>::max_num_args_called; 
  static constexpr bool is_timedep = traits<XprType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<XprType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = traits<XprType>::maxOrder; // highest order of derivative in the expression 
  typedef typename traits<XprType>::orders orders;
  typedef typename traits<XprType>::node_selector_tag node_selector_tag; // give priority to lhs 
}; 

// map to base tag 
template<class UnaryOp, class XprType>
struct map_to_base_tag<fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>> : map_to_base_tag<std::remove_cv_t<std::remove_reference_t<XprType>>>{};

} // end namespace internal 
} // end namespace linops
} // end namespace fornfdm 

namespace Eigen{
namespace internal{

template<class UnaryOp, class XprType>
struct traits<fornfdm::linops::NwiseUnaryOp<UnaryOp,XprType>>
{
  typedef fornfdm::Scalar Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Sparse StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = Eigen::RowMajorBit | NestByRefBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<class UnaryOp, class _XprType>
struct evaluator<fornfdm::linops::NwiseUnaryOp<UnaryOp,_XprType>> 
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::NwiseUnaryOp<UnaryOp,_XprType>>
{
  using XprType = fornfdm::linops::NwiseUnaryOp<UnaryOp,_XprType>; 
  using Impl = typename fornfdm::linops::internal::EigenEvaluatorImpl<XprType>; 
  enum {CoeffReadCost = Impl::CoeffReadCost, Flags = Impl::Flags};
  using InnerIterator = typename Impl::InnerIterator; 
  evaluator(const XprType& xpr)
    : Impl(xpr)
  {}
}; 

} // end namespac internal 
} // end namespac Eigen

namespace fornfdm{ 
namespace linops{

template<class UnaryOp, class XprType>
class NwiseUnaryOp : public fornfdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fornfdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(NwiseUnaryOp)
    typedef typename fornfdm::linops::internal::NestedStorage<XprType>::type XprTypeNested;
    typedef typename std::remove_reference<std::remove_cv_t<XprType>>::type NestedExpression;

  protected:
    // Member data ----------------------------------- 
    UnaryOp m_functor; 
    XprTypeNested m_xpr; 
  
  public:
    // Constructors ====================== 
    NwiseUnaryOp(XprType&& xpr, UnaryOp func = UnaryOp())
      : m_xpr(std::forward<XprType>(xpr)), m_functor(func) 
    {}

    // Member Funcs ----------------------------------- 
    const auto& functor() const { return m_functor; }
    const auto& nestedExpression() const { return m_xpr; }
    auto& nestedExpression(){ return m_xpr; }
}; 

template<typename XprType, typename = std::enable_if_t<fornfdm::linops::internal::is_partialderiv_crtp<XprType>::value> >
auto operator-(XprType&& xpr)
{
  return NwiseUnaryOp<fornfdm::linops::internal::UnaryNegateFO, XprType>(std::forward<XprType>(xpr), fornfdm::linops::internal::UnaryNegateFO{}); 
}; 

template<typename C, typename XprType, typename = std::enable_if_t<std::is_convertible<C,fornfdm::Scalar>::value && fornfdm::linops::internal::is_partialderiv_crtp<XprType>::value>>
auto operator*(C&& c, XprType&& xpr)
{
  return NwiseUnaryOp<fornfdm::linops::internal::UnaryScalarMultiplyFO<C>, XprType>(std::forward<XprType>(xpr), fornfdm::linops::internal::UnaryScalarMultiplyFO<C>(std::forward<C>(c))); 
}; 

} // end namespace linops 
} // end namespace fornfdm 

#endif // NwiseUnaryOp.hpp