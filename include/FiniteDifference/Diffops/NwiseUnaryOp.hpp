// NwiseUnaryOp.hpp
//
//
//
// JAF 4/22/2026 

#ifndef DIFFOPS_NWISEUNARYOP_H
#define DIFFOPS_NWISEUNARYOP_H

#include "EvaluatorBase.hpp"
#include "Functors.hpp"

namespace fdm{
namespace linops{

// Forward Declaration
template<class UnaryOp, class XprType>
struct NwiseUnaryOp; 

namespace internal{

// Evaluator 
template<class UnaryOp, class XprType>
struct Evaluator<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>> : public EvaluatorBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>
{
  const fdm::linops::NwiseUnaryOp<UnaryOp,XprType>& m_xpr; 
  Evaluator<XprType> m_nested_eval; 
  Evaluator(const fdm::linops::NwiseUnaryOp<UnaryOp,XprType>& xpr) : m_xpr(xpr), m_nested_eval(xpr.nestedExpression()){} 
  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const Coordinate<N>& coords) const 
  {
    return m_xpr.functor()( m_nested_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords)); 
  }
};

// Traits
template<class UnaryOp, class XprType>
struct traits_impl<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>> : public traits_impl<XprType>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = true; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
}; 

} // end namespace internal 
} // end namespace linops
} // end namespace fdm 

namespace Eigen{
namespace internal{

template<class UnaryOp, class XprType>
struct traits<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>
{
  typedef fdm::Scalar Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Sparse StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = Eigen::RowMajor,  /* no | NestByRefBit */ /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

} // end namespac internal 
} // end namespac internal 

namespace fdm{ 
namespace linops{

template<class UnaryOp, class XprType>
class NwiseUnaryOp : public fdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(NwiseUnaryOp)
    typedef typename Eigen::internal::ref_selector<XprType>::type XprTypeNested;
    typedef typename std::remove_reference<std::remove_cv_t<XprType>>::type NestedExpression;

  protected:
    // Member data ----------------------------------- 
    UnaryOp m_functor; 
    XprTypeNested m_xpr; 
  
  public:
    // Constructors ====================== 
    NwiseUnaryOp(const XprType& xpr, const UnaryOp& func = UnaryOp())
      : m_xpr(xpr), m_functor(func) 
    {}

    // Member Funcs ----------------------------------- 
    const auto& functor() const { return m_functor; }
    const auto& nestedExpression() const { return m_xpr; }
    auto& nestedExpression(){ return m_xpr; }
}; 

template<typename Derived>
auto operator-(const PartialDerivBase<Derived>& xpr)
{
  return NwiseUnaryOp(xpr.derived(), fdm::linops::internal::UnaryNegateFO{}); 
}; 

template<typename C, typename Derived>
auto operator*(C&& c, const PartialDerivBase<Derived>& xpr)
{
  return NwiseUnaryOp(xpr.derived(), fdm::linops::internal::UnaryScalarMultiplyFO<C>(std::forward<C>(c))); 
}; 

} // end namespace linops 
} // end namespace fdm 

#endif 