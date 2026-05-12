// NwiseBinaryOp.hpp
//
//
//
// JAF 4/23/2026 

#ifndef FDM_DIFFOPS_NWISEBINARYOP_H
#define FDM_DIFFOPS_NWISEBINARYOP_H

#include "../Traits.hpp"
#include "Traits.hpp"
#include "EvaluatorBase.hpp"
#include "Functors.hpp"

namespace fdm{
namespace linops{

// Forward Declaration
template<class BinaryOp, class LhsType, class RhsType>
struct NwiseBinaryOp; 

namespace internal{ 

// Evaluator 
template<class BinaryOp, class LhsType, class RhsType>
struct Evaluator<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>> : public EvaluatorBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  using XprType = fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>; 
  const XprType& m_xpr; 
  Evaluator<typename XprType::Lhs> m_lhs_eval; 
  Evaluator<typename XprType::Rhs> m_rhs_eval; 
  Evaluator(const fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>& xpr)
    : m_xpr(xpr), m_lhs_eval(xpr.lhs()), m_rhs_eval(xpr.rhs())
  {}

  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const fdm::Coordinate<N>& coords) const 
  {
    return m_xpr.functor()( 
      m_lhs_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords), 
      m_rhs_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords)
    ); 
  }
};

// Traits
template<class BinaryOp, class LhsType, class RhsType>
struct traits_impl<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = std::max(traits<LhsType>::max_num_args_called,traits<RhsType>::max_num_args_called); 
  static constexpr bool is_timedep = traits<LhsType>::is_timedep || traits<RhsType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<LhsType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t maxOrder = std::max(traits<LhsType>::maxOrder,traits<RhsType>::maxOrder); // highest order of derivative in the expression 
  typedef typename promote_node_selector_type<typename traits<LhsType>::node_selector_tag,typename traits<RhsType>::node_selector_tag>::type node_selector_tag; // gurantees both minimum are fulfilled
}; 

} // end namespace internal 
} // end namespace linops
} // end namespace fdm 

namespace Eigen{
namespace internal{

template<class BinaryOp, class LhsType, class RhsType>
struct traits<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
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
    Flags = Eigen::RowMajorBit | NestByRefBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<class BinaryOp, class LhsType, class RhsType>
struct evaluator<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>> 
  : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>, 
  public evaluator_base<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>> 
{
  struct InnerIterator
    : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>::InnerIterator
  {
    enum { 
      CoeffReadCost = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>::CoeffReadCost, 
      Flags = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>::Flags 
    };

    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>::InnerIterator(eval, row_idx)
    {}
  }; 
  evaluator(const fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>& xpr_d)
    : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NwiseBinaryOp<BinaryOp,LhsType,RhsType>>>(xpr_d)
  {}
}; 

} // end namespac internal 
} // end namespac Eigen 

namespace fdm{ 
namespace linops{

template<class BinaryOp, class LhsType, class RhsType>
class NwiseBinaryOp : public fdm::linops::PartialDerivBase<NwiseBinaryOp<BinaryOp,LhsType,RhsType>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fdm::linops::PartialDerivBase<NwiseBinaryOp<BinaryOp,LhsType,RhsType>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(NwiseBinaryOp)
    typedef typename std::remove_cv_t<std::remove_reference_t<LhsType>> Lhs; 
    typedef typename std::remove_cv_t<std::remove_reference_t<RhsType>> Rhs; 
    typedef typename fdm::linops::internal::NestedStorage<LhsType>::type LhsNested;
    typedef typename fdm::linops::internal::NestedStorage<RhsType>::type RhsNested;

    // Friends 
    friend Eigen::internal::evaluator<NwiseBinaryOp>; 
    friend fdm::linops::internal::EvaluatorBase<NwiseBinaryOp>; 
    friend fdm::linops::internal::Evaluator<NwiseBinaryOp>;

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
    void setTime_hooked(double t)
    {
      m_lhs.const_cast_derived().setTime_hooked(t); 
      m_rhs.const_cast_derived().setTime_hooked(t); 
    }
}; 

template<
  typename LeftArg,
  typename RightArg, 
  typename = std::enable_if_t<
    fdm::linops::internal::is_partialderiv_crtp<LeftArg>::value &&
    fdm::linops::internal::is_partialderiv_crtp<RightArg>::value &&  
    fdm::internal::is_matching<
      typename linops::internal::traits<LeftArg>::node_selector_tag, 
      typename linops::internal::traits<RightArg>::node_selector_tag
    >::value &&  
    (linops::internal::traits<LeftArg>::is_timedep == linops::internal::traits<RightArg>::is_timedep) && 
    (linops::internal::traits<LeftArg>::direction == linops::internal::traits<RightArg>::direction)
  >
>
auto operator-(LeftArg&& lhs, RightArg&& rhs)
{
  return NwiseBinaryOp<fdm::linops::internal::BinarySubtractionFO, LeftArg, RightArg>(std::forward<LeftArg>(lhs), std::forward<RightArg>(rhs), fdm::linops::internal::BinarySubtractionFO{}); 
}; 

template<
  typename LeftArg,
  typename RightArg, 
  typename = std::enable_if_t<
    fdm::linops::internal::is_partialderiv_crtp<LeftArg>::value &&
    fdm::linops::internal::is_partialderiv_crtp<RightArg>::value &&  
    fdm::internal::is_matching<
      typename linops::internal::traits<LeftArg>::node_selector_tag, 
      typename linops::internal::traits<RightArg>::node_selector_tag
    >::value &&  
    (linops::internal::traits<LeftArg>::is_timedep == linops::internal::traits<RightArg>::is_timedep) && 
    (linops::internal::traits<LeftArg>::direction == linops::internal::traits<RightArg>::direction)
  >
>
auto operator+(LeftArg&& lhs, RightArg&& rhs)
{
  return NwiseBinaryOp<fdm::linops::internal::BinaryAdditionFO, LeftArg, RightArg>(std::forward<LeftArg>(lhs), std::forward<RightArg>(rhs), fdm::linops::internal::BinaryAdditionFO{}); 
}; 

} // end namespace linops 
} // end namespace fdm 

#endif // NwiseBinaryOp.hpp