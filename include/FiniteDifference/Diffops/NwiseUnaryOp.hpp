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
  using UnarXprType = fdm::linops::NwiseUnaryOp<UnaryOp,XprType>;  
  const UnarXprType& m_xpr; 
  Evaluator<typename UnarXprType::NestedExpression> m_nested_eval; 
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
  static constexpr std::size_t max_num_args_called = traits<XprType>::max_num_args_called; 
  static constexpr bool is_timedep = traits<XprType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<XprType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = traits<XprType>::maxOrder; // highest order of derivative in the expression 
  typedef typename traits<XprType>::node_selector_tag node_selector_tag; // give priority to lhs 
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
    Flags = Eigen::RowMajorBit,  /* no | NestByRefBit */ /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<class UnaryOp, class XprType>
struct evaluator<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>> 
  : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>, 
  public evaluator_base<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>> 
{
  struct InnerIterator
    : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>::InnerIterator
  {
    enum { 
      CoeffReadCost = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>::CoeffReadCost, 
      Flags = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>::Flags 
    };

    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>::InnerIterator(eval, row_idx)
    {}
  }; 
  evaluator(const fdm::linops::NwiseUnaryOp<UnaryOp,XprType>& xpr_d)
    : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseUnaryOp<UnaryOp,XprType>>>(xpr_d)
  {}
}; 

} // end namespac internal 
} // end namespac Eigen

namespace fdm{ 
namespace linops{

template<class UnaryOp, class XprType>
class NwiseUnaryOp : public fdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fdm::linops::PartialDerivBase<NwiseUnaryOp<UnaryOp,XprType>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(NwiseUnaryOp)
    typedef typename fdm::linops::internal::NestedStorage<XprType>::type XprTypeNested;
    typedef typename std::remove_reference<std::remove_cv_t<XprType>>::type NestedExpression;

    // Friends ------------------------------- 
    friend Eigen::internal::evaluator<NwiseUnaryOp>; 
    friend fdm::linops::internal::EvaluatorBase<NwiseUnaryOp>; 
    friend fdm::linops::internal::Evaluator<NwiseUnaryOp>;

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
    void setTime_hooked(double t)
    {
      m_xpr.const_cast_derived().setTime_hooked(t); 
    }
}; 

template<typename XprType, typename = std::enable_if_t<fdm::linops::internal::is_partialderiv_crtp<XprType>::value> >
auto operator-(XprType&& xpr)
{
  return NwiseUnaryOp<fdm::linops::internal::UnaryNegateFO, XprType>(std::forward<XprType>(xpr), fdm::linops::internal::UnaryNegateFO{}); 
}; 

template<typename C, typename XprType, typename = std::enable_if_t<std::is_arithmetic<C>::value && fdm::linops::internal::is_partialderiv_crtp<XprType>::value>>
auto operator*(C&& c, XprType&& xpr)
{
  return NwiseUnaryOp<fdm::linops::internal::UnaryScalarMultiplyFO<C>, XprType>(std::forward<XprType>(xpr), fdm::linops::internal::UnaryScalarMultiplyFO<C>(std::forward<C>(c))); 
}; 

} // end namespace linops 
} // end namespace fdm 

#endif 