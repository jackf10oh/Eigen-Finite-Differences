// NwiseBinaryOp.hpp
//
//
//
// JAF 4/23/2026 

#ifndef FORNFDM_DIFFOPS_NWISEBINARYOP_H
#define FORNFDM_DIFFOPS_NWISEBINARYOP_H

#include "../types.hpp"
#include "../traits.hpp"
#include "traits.hpp"
#include "EvaluatorBase.hpp"
#include "PartialDerivBase.hpp"
#include "EigenEvaluatorImpl.hpp"
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
  Evaluator(const fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>& xpr, fornfdm::Real t)
    : EvaluatorBase<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>(t), 
    m_xpr(xpr), 
    m_lhs_eval(xpr.lhs(),t), 
    m_rhs_eval(xpr.rhs(),t)
  {}

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return [f = m_xpr.functor(), n1 = m_lhs_eval.createReader(coord,t), n2 = m_rhs_eval.createReader(coord,t)](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride)
    {
      return f(n1(weights,idx,stride), n2(weights,idx,stride));
    };
  }

  template<std::size_t N>
  struct ExactReader
  {
    using LeftNestedReader = decltype(std::declval<const Evaluator<typename XprType::Lhs>&>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>(), std::declval<fornfdm::Real>()));
    using RightNestedReader = decltype(std::declval<const Evaluator<typename XprType::Rhs>&>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>(), std::declval<fornfdm::Real>()));
    const XprType& m_xpr;
    LeftNestedReader m_nested_left;
    RightNestedReader m_nested_right;

    ExactReader(const Evaluator& eval, const fornfdm::Coordinate<N>& coord, fornfdm::Real t)
      : m_xpr(eval.m_xpr), 
      m_nested_left(eval.m_lhs_eval.template createExactReader<N>(coord,t)),
      m_nested_right(eval.m_rhs_eval.template createExactReader<N>(coord,t))
    {}

    template<std::size_t... orders>
    fornfdm::Scalar operator()(const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride, std::index_sequence<orders...>) const
    {
      return m_xpr.functor()(m_nested_left(weights, idx, stride, std::index_sequence<orders...>{}), m_nested_right(weights, idx, stride, std::index_sequence<orders...>{}));
    }
  };

  template<std::size_t N>
  auto createExactReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return ExactReader<N>(*this, coord, t);
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
  static constexpr std::size_t max_arity = std::max(traits<LhsType>::max_arity,traits<RhsType>::max_arity); 
  static constexpr bool is_timedep = traits<LhsType>::is_timedep || traits<RhsType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<LhsType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t max_order = std::max(traits<LhsType>::max_order,traits<RhsType>::max_order); // highest order of derivative in the expression 
  typedef typename promote_node_selector_tags<typename traits<LhsType>::node_selector_tag,typename traits<RhsType>::node_selector_tag>::type node_selector_tag; // gurantees both minimum are fulfilled
  typedef typename merge_orders< typename traits<LhsType>::orders, typename traits<RhsType>::orders>::type orders;
};

// ==============================
// map_to_base specializations 
// ==============================

// (direction 0 && max_num_args <= 1) --> LeftKroneckerTag or TimeDepLeftKroneckerTag
template<class BinaryOp, class LhsType, class RhsType>
struct map_to_base_tag<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>,
  std::enable_if_t<
    traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::direction == 0 && 
    traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::max_arity <= 1 
  >
>
{
  using traits_t = traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>;
  using type = typename std::conditional<traits_t::is_timedep, TimeDepLeftKroneckerTag, LeftKroneckerTag>::type;
};

// (direction 0 && max_num_args > 1) --> StoredWeightsTag or TimeDepStoredWeightsTag
template<class BinaryOp, class LhsType, class RhsType>
struct map_to_base_tag<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>,
  std::enable_if_t<
    (traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::direction == 0 && 
    traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::max_arity > 1)
  >
>
{
  using traits_t = traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>;
  #ifndef FORNFDM_STORE_FULL_KRONECKER
  using type = typename std::conditional<traits_t::is_timedep, TimeDepStoredWeightsTag, StoredWeightsTag>::type;
  #else
  // Goes to a LeftKronecker instead!
  using type = typename std::conditional<traits_t::is_timedep, TimeDepLeftKroneckerTag, LeftKroneckerTag>::type;
  #endif
};

// (direction != 0 && max_num_args == 0) --> DoubleKroneckerTag or TimeDepDoubleKroneckerTag
template<class BinaryOp, class LhsType, class RhsType>
struct map_to_base_tag<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>,
  std::enable_if_t<
    (traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::direction != 0 && 
    traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::max_arity == 0)
  >
>
{
  using traits_t = traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>;
  using type = typename std::conditional<traits_t::is_timedep, TimeDepDoubleKroneckerTag, DoubleKroneckerTag>::type;
};

// (direction != 0 && max_num_args != 0) --> StoredWeightsTag or TimeDepStoredWeightsTag
template<class BinaryOp, class LhsType, class RhsType>
struct map_to_base_tag<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>,
  std::enable_if_t<
    (traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::direction != 0 && 
    traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>::max_arity != 0)
  >
>
{
  using traits_t = traits<fornfdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>;
  #ifndef FORNFDM_STORE_FULL_KRONECKER
  using type = typename std::conditional<traits_t::is_timedep, TimeDepStoredWeightsTag, StoredWeightsTag>::type;
  #else
  // Goes to a LeftKronecker instead!
  using type = typename std::conditional<traits_t::is_timedep, TimeDepLeftKroneckerTag, LeftKroneckerTag>::type;
  #endif
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
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>>
{
  // Flags ---
  using XprType = fornfdm::linops::NwiseBinaryOp<BinaryOp, LhsType, RhsType>; 
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