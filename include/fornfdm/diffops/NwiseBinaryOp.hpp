// NwiseBinaryOp.hpp
//
//
//
// JAF 4/23/2026 

#ifndef FORNFDM_DIFFOPS_NWISEBINARYOP_H
#define FORNFDM_DIFFOPS_NWISEBINARYOP_H

#include "../traits.hpp"
#include "traits.hpp"
#include "KroneckerEvaluator.hpp"
#include "EvaluatorBase.hpp"
#include "functors.hpp"

namespace fornfdm{
namespace linops{

// Forward Declaration
template<class BinaryOp, class LhsType, class RhsType>
struct NwiseBinaryOp; 

namespace internal{ 

// Evaluator 
template<class BinaryOp, class LhsType, class RhsType>
struct Evaluator<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>> : public EvaluatorBase<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  using XprType = fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>; 
  const XprType& m_xpr; 
  Evaluator<typename XprType::Lhs> m_lhs_eval; 
  Evaluator<typename XprType::Rhs> m_rhs_eval; 
  Evaluator(const fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>& xpr)
    : m_xpr(xpr), m_lhs_eval(xpr.lhs()), m_rhs_eval(xpr.rhs())
  {}

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return [f = m_xpr.functor(), n1 = m_lhs_eval.createReader(coord,t), n2 = m_rhs_eval.createReader(coord,t)](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride)
    {
      return f(n1(weights,idx,stride), n2(weights,idx,stride));
    };
  }
};

// traits
template<class BinaryOp, class LhsType, class RhsType>
struct traits_impl<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = std::max(traits<LhsType>::max_num_args_called,traits<RhsType>::max_num_args_called); 
  static constexpr bool is_timedep = traits<LhsType>::is_timedep || traits<RhsType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<LhsType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = std::max(traits<LhsType>::maxOrder,traits<RhsType>::maxOrder); // highest order of derivative in the expression 
  typedef typename promote_node_selector_tags<typename traits<LhsType>::node_selector_tag,typename traits<RhsType>::node_selector_tag>::type node_selector_tag; // gurantees both minimum are fulfilled
}; 

} // end namespace internal 
} // end namespace linops
} // end namespace fornfdm 

namespace Eigen{
namespace internal{

template<class BinaryOp, class LhsType, class RhsType>
struct traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
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

template<class BinaryOp, class LhsType, class RhsType>
struct evaluator<fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>> 
  : public evaluator_base<fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>>, 
  public fornfdm::linops::internal::KroneckerEvaluator<fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>>
{
  using XprType = fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>; 
  using Impl = typename fornfdm::linops::internal::KroneckerEvaluator<fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>>; 
  using InnerIterator = typename Impl::InnerIterator; 
  evaluator(const XprType& xpr)
    : Impl(xpr)
  {}
};

} // end namespac internal 
} // end namespac Eigen 

namespace fornfdm{ 
namespace linops{

template<class BinaryOp, class LhsType, class RhsType>
class NwiseBinaryOp : public fornfdm::linops::PartialDerivBase<NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fornfdm::linops::PartialDerivBase<NwiseBinaryOp<BinaryOp,LhsType,RhsType>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(NwiseBinaryOp)
    typedef typename std::remove_cv_t<std::remove_reference_t<LhsType>> Lhs; 
    typedef typename std::remove_cv_t<std::remove_reference_t<RhsType>> Rhs; 
    typedef typename fornfdm::linops::internal::NestedStorage<LhsType>::type LhsNested;
    typedef typename fornfdm::linops::internal::NestedStorage<RhsType>::type RhsNested;

    // Friends 
    friend Eigen::internal::evaluator<NwiseBinaryOp>; 
    friend fornfdm::linops::internal::KroneckerEvaluator<NwiseBinaryOp>; 
    friend fornfdm::linops::internal::EvaluatorBase<NwiseBinaryOp>; 
    friend fornfdm::linops::internal::Evaluator<NwiseBinaryOp>;

  protected:
    // Member data ----------------------------------- 
    BinaryOp m_functor; 
    LhsNested m_lhs; 
    RhsNested m_rhs; 
  
  public:
    // Constructors ====================== 
    NwiseBinaryOp(LhsType&& lhs, RhsType&& rhs, BinaryOp func = BinaryOp())
      : m_lhs(std::forward<LhsType>(lhs)), m_rhs(std::forward<RhsType>(rhs)), m_functor(func) 
    {}

    // Member Funcs ----------------------------------- 
    const auto& functor() const { return m_functor; }
    const auto& lhs() const { return m_lhs; }
    auto& lhs(){ return m_lhs; }
    const auto& rhs() const { return m_rhs; }
    auto& rhs(){ return m_rhs; }
}; 

template<
  typename LeftArg,
  typename RightArg, 
  typename = std::enable_if_t<
    fornfdm::linops::internal::is_partialderiv_crtp<LeftArg>::value &&
    fornfdm::linops::internal::is_partialderiv_crtp<RightArg>::value &&  
    (linops::internal::traits<LeftArg>::direction == linops::internal::traits<RightArg>::direction) &&
    fornfdm::linops::internal::promote_node_selector_tags<
      typename linops::internal::traits<LeftArg>::node_selector_tag, 
      typename linops::internal::traits<RightArg>::node_selector_tag
    >::is_match
  >
>
auto operator-(LeftArg&& lhs, RightArg&& rhs)
{
  return NwiseBinaryOp<fornfdm::linops::internal::BinarySubtractionFO, LeftArg, RightArg>(std::forward<LeftArg>(lhs), std::forward<RightArg>(rhs), fornfdm::linops::internal::BinarySubtractionFO{}); 
}; 

template<
  typename LeftArg,
  typename RightArg, 
  typename = std::enable_if_t<
    internal::is_partialderiv_crtp<LeftArg>::value &&
    internal::is_partialderiv_crtp<RightArg>::value &&
    (linops::internal::traits<LeftArg>::direction == linops::internal::traits<RightArg>::direction) &&
    internal::promote_node_selector_tags<
      typename internal::traits<LeftArg>::node_selector_tag,
      typename internal::traits<RightArg>::node_selector_tag
    >::is_match
  >
>
auto operator+(LeftArg&& lhs, RightArg&& rhs)
{
  return NwiseBinaryOp<fornfdm::linops::internal::BinaryAdditionFO, LeftArg, RightArg>(std::forward<LeftArg>(lhs), std::forward<RightArg>(rhs), fornfdm::linops::internal::BinaryAdditionFO{}); 
}; 

} // end namespace linops 
} // end namespace fornfdm 

#endif // NwiseBinaryOp.hpp