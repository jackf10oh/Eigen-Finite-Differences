// Product.hpp 
//
// Binary Expression of CoeffBase<> * PartialDerivBase<> 
//
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_COEFFPRODUCT_H
#define FORNFDM_COEFFS_COEFFPRODUCT_H

#include "../types.hpp"
#include "../diffops/traits.hpp"
#include "../diffops/functors.hpp"
#include "../diffops/EvaluatorBase.hpp"
#include "../diffops/PartialDerivBase.hpp"
#include "../diffops/EigenEvaluatorImpl.hpp"

namespace fornfdm{
namespace linops{

// forward declaration 
template<class LeftCoeff, class RightDeriv>
class CoeffProduct; 

namespace internal{

// traits 
template<class LeftCoeff, class RightDeriv>
struct traits_impl<CoeffProduct<LeftCoeff, RightDeriv>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = true; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = std::max(traits<LeftCoeff>::max_num_args_called,traits<RightDeriv>::max_num_args_called); 
  static constexpr bool is_timedep = traits<LeftCoeff>::is_timedep || traits<RightDeriv>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<RightDeriv>::direction; // give priority to RHS 
  static constexpr std::size_t maxOrder = traits<RightDeriv>::maxOrder; // highest order of derivative in the expression 
  typedef typename traits<RightDeriv>::node_selector_tag node_selector_tag; // give priority to RHS 
  typedef typename traits<RightDeriv>::orders orders; // LeftCoeff never has any orders 
}; 

// Evaluator
template<class LeftCoeff, class RightDeriv>
struct Evaluator<CoeffProduct<LeftCoeff, RightDeriv>> : public EvaluatorBase<CoeffProduct<LeftCoeff, RightDeriv>>
{
  using XprType = fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>; 
  const XprType& m_xpr; 
  Evaluator<typename XprType::Rhs> m_rhs_eval; 
  Evaluator(const fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>& xpr)
    : m_xpr(xpr), m_rhs_eval(xpr.rhs())
  {}

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    if constexpr(fornfdm::linops::internal::traits<LeftCoeff>::is_timedep){
      return [c = coord.applyBindFirst(m_xpr.lhs().callable(), t), nested = m_rhs_eval.createReader(coord,t)](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride)
      {
        return c * nested(weights, idx, stride);
      };
    }
    else{
      return [c = coord.apply(m_xpr.lhs().callable()), nested = m_rhs_eval.createReader(coord,t)](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride)
      {
        return c * nested(weights, idx, stride);
      };
    }
  }

  template<std::size_t N>
  struct AutonomousExactReader
  {
    using RhsReader = decltype(std::declval<const Evaluator<typename XprType::Rhs>&>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>(), std::declval<fornfdm::Real>()));    
    AutonomousExactReader(const Evaluator& eval, const fornfdm::Coordinate<N>& coord, fornfdm::Real t)
      : m_rhs_reader(eval.m_rhs_eval.template createExactReader<N>(coord, t)),
      m_coeff(coord.apply(eval.m_xpr.lhs().callable()))
    {}
    RhsReader m_rhs_reader;
    fornfdm::Scalar m_coeff; 
    template<std::size_t... orders>
    fornfdm::Scalar operator()(const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride, std::index_sequence<orders...>) const
    {
      return m_coeff * m_rhs_reader(weights, idx, stride, std::index_sequence<orders...>{});
    }
  };

  template<std::size_t N>
  struct TimeDepExactReader
  {
    using RhsReader = decltype(std::declval<const Evaluator<typename XprType::Rhs>&>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>(), std::declval<fornfdm::Real>()));    
    TimeDepExactReader(const Evaluator& eval, const fornfdm::Coordinate<N>& coord, fornfdm::Real t)
      : m_rhs_reader(eval.m_rhs_eval.template createExactReader<N>(coord, t)),
      m_coeff(coord.applyBindFirst(eval.m_xpr.lhs().callable(), t))
    {}
    RhsReader m_rhs_reader;
    fornfdm::Scalar m_coeff; 
    template<std::size_t... orders>
    fornfdm::Scalar operator()(const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride, std::index_sequence<orders...>) const
    {
      return m_coeff * m_rhs_reader(weights, idx, stride, std::index_sequence<orders...>{});
    }
  };

  template<std::size_t N>
  auto createExactReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    if constexpr(traits<LeftCoeff>::is_timedep){
      return TimeDepExactReader<N>(*this, coord, t);
    }
    else{
      return AutonomousExactReader<N>(*this, coord, t);
    }
  }
}; 

// ==============================
// map_to_base specializations 
// ==============================

// (direction 0 && max_num_args <= 1) --> LeftKroneckerTag or TimeDepLeftKroneckerTag
template<class LeftCoeff, class RightDeriv>
struct map_to_base_tag<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>,
  std::enable_if_t<
    (traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::direction == 0 && 
    traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::max_num_args_called <= 1)
  >
>
{
  using traits_t = traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>;
  using type = typename std::conditional<traits_t::is_timedep, TimeDepLeftKroneckerTag, LeftKroneckerTag>::type;
};

// (direction != 0 && max_num_args == 0) --> DoubleKroneckerTag or TimeDepDoubleKroneckerTag
template<class LeftCoeff, class RightDeriv>
struct map_to_base_tag<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>,
  std::enable_if_t<
    (traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::direction != 0 && 
    traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::max_num_args_called == 0)
  >
>
{
  using traits_t = traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>;
  using type = typename std::conditional<traits_t::is_timedep, TimeDepDoubleKroneckerTag, DoubleKroneckerTag>::type;
};

// (direction 0 && max_num_args > 1) --> StoredWeightsTag or TimeDepStoredWeightsTag
template<class LeftCoeff, class RightDeriv>
struct map_to_base_tag<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>,
  std::enable_if_t<
    (traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::direction == 0 && 
    traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::max_num_args_called > 1)
  >
>
{
  using traits_t = traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>;
  #ifndef FORNFDM_STORE_FULL_KRONECKER
  using type = typename std::conditional<traits_t::is_timedep, TimeDepStoredWeightsTag, StoredWeightsTag>::type;
  #else
  // Goes to a LeftKronecker instead!
  using type = typename std::conditional<traits_t::is_timedep, TimeDepLeftKroneckerTag, LeftKroneckerTag>::type;
  #endif
};

// (direction != 0 && max_num_args != 0) --> StoredWeightsTag or TimeDepStoredWeightsTag
template<class LeftCoeff, class RightDeriv>
struct map_to_base_tag<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>,
  std::enable_if_t<
    (traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::direction != 0 && 
    traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>::max_num_args_called != 0)
  >
>
{
  using traits_t = traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>;
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

// traits
template<class LeftCoeff, class RightDeriv>
struct traits<fornfdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>
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
    Flags = Eigen::RowMajorBit| NestByRefBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<class LeftCoeff, class RightDeriv>
struct evaluator<fornfdm::linops::CoeffProduct<LeftCoeff,RightDeriv>> 
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>
{
  using XprType = fornfdm::linops::CoeffProduct<LeftCoeff,RightDeriv>; 
  using Impl = typename fornfdm::linops::internal::EigenEvaluatorImpl<XprType>; 
  using InnerIterator = typename Impl::InnerIterator;
  evaluator(const XprType& xpr)
    : Impl(xpr)
  {}
};

} // end namespace internal 
} // end namespace Eigen 

namespace fornfdm{
namespace linops{ 

template<class LeftCoeff, class RightDeriv>
class CoeffProduct : public fornfdm::linops::PartialDerivBase<CoeffProduct<LeftCoeff, RightDeriv>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fornfdm::linops::PartialDerivBase<CoeffProduct<LeftCoeff, RightDeriv>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(CoeffProduct)
    typedef typename std::remove_cv_t<std::remove_reference_t<LeftCoeff>> Lhs; 
    typedef typename std::remove_cv_t<std::remove_reference_t<RightDeriv>> Rhs; 
    typedef typename fornfdm::linops::internal::NestedStorage<LeftCoeff>::type LhsNested;
    typedef typename fornfdm::linops::internal::NestedStorage<RightDeriv>::type RhsNested;

  protected:
    // Member data ----------------------------------- 
    fornfdm::linops::internal::BinaryMultiplyFO m_functor; 
    LhsNested m_lhs; 
    RhsNested m_rhs; 
  
  public:
    // Constructors ====================== 
    CoeffProduct(LeftCoeff&& lhs, RightDeriv&& rhs)
      : m_lhs(lhs), m_rhs(rhs), m_functor() 
    {}

    // Member Funcs ----------------------------------- 
    const auto& functor() const { return m_functor; }
    const auto& lhs() const { return m_lhs; }
    auto& lhs(){ return m_lhs; }
    const auto& rhs() const { return m_rhs; }
    auto& rhs(){ return m_rhs; }
}; 

} // end namespace linops  
} // end namespace fornfdm 

#endif // CoeffProduct.hpp 