// Product.hpp 
//
// Binary Expression of CoeffBase<> * PartialDerivBase<> 
//
// JAF 4/24/2026 

#ifndef FDM_COEFF_PRODUCT_H
#define FDM_COEFF_PRODUCT_H

#include "../Diffops/EvaluatorBase.hpp"

namespace fdm{
namespace linops{

// forward declaration 
template<class LeftCoeff, class RightDeriv>
class CoeffProduct; 

namespace internal{

// Traits 
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
}; 

// Evaluator
template<class LeftCoeff, class RightDeriv>
struct Evaluator<CoeffProduct<LeftCoeff, RightDeriv>> : public EvaluatorBase<CoeffProduct<LeftCoeff, RightDeriv>>
{
  const fdm::linops::CoeffProduct<LeftCoeff, RightDeriv>& m_xpr; 
  Evaluator<RightDeriv> m_rhs_eval; 
  Evaluator(const fdm::linops::CoeffProduct<LeftCoeff, RightDeriv>& xpr)
    : m_xpr(xpr), m_rhs_eval(xpr.rhs())
  {}

  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const Coordinate<N>& coords) const 
  {
    return m_xpr.functor()(coords.apply(m_xpr.lhs().callable()), m_rhs_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords)); 
  }
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm 

namespace Eigen{
namespace internal{

// Traits
template<class LeftCoeff, class RightDeriv>
struct traits<fdm::linops::CoeffProduct<LeftCoeff, RightDeriv>>
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

// Evaluator 
template<class LeftCoeff, class RightDeriv>
struct evaluator<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>> 
  : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>, 
  public evaluator_base<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>> 
{
  struct InnerIterator
    : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>::InnerIterator
  {
    enum { 
      CoeffReadCost = evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>::CoeffReadCost, 
      Flags = evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>::Flags 
    };

    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>::InnerIterator(eval, row_idx)
    {}
  }; 
  evaluator(const fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>& xpr_d)
    : evaluator<fdm::linops::PartialDerivBase<fdm::linops::CoeffProduct<LeftCoeff,RightDeriv>>>(xpr_d)
  {}
}; 

} // end namespace internal 
} // end namespace Eigen 

namespace fdm{
namespace linops{ 

template<class LeftCoeff, class RightDeriv>
class CoeffProduct : public fdm::linops::PartialDerivBase<CoeffProduct<LeftCoeff, RightDeriv>>
{
  public: 
    // Type Defs ------------------------------------- 
    using Base = fdm::linops::PartialDerivBase<CoeffProduct<LeftCoeff, RightDeriv>>; 
    EIGEN_SPARSE_PUBLIC_INTERFACE(CoeffProduct)
    typedef typename std::remove_cv_t<std::remove_reference_t<LeftCoeff>> Lhs; 
    typedef typename std::remove_cv_t<std::remove_reference_t<RightDeriv>> Rhs; 
    typedef typename Eigen::internal::ref_selector<LeftCoeff>::type LhsNested;
    typedef typename Eigen::internal::ref_selector<RightDeriv>::type RhsNested;

  protected:
    // Member data ----------------------------------- 
    fdm::linops::internal::BinaryMultiplyFO m_functor; 
    LhsNested m_lhs; 
    RhsNested m_rhs; 
  
  public:
    // Constructors ====================== 
    CoeffProduct(const LeftCoeff& lhs, const RightDeriv& rhs)
      : m_lhs(lhs), m_rhs(rhs), m_functor() 
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

} // end namespace linops  
} // end namespace fdm 

#endif 