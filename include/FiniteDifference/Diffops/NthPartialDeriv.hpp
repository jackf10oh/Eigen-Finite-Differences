// NthPartialDeriv.hpp 
//
// Concrete class for Nth Partial derivative in direction. 
// 
// JAF 4/17/2026 

#ifndef FDM_DIFFOPS_NTHPARTIALDERIV_H
#define FDM_DIFFOPS_NTHPARTIALDERIV_H

#include "EvaluatorBase.hpp"
#include "CenteredNodeSelector.hpp"
#include "../Traits.hpp"

namespace fdm{
namespace linops{

  // Forward Declaration ------------------------------------------------- 
template<std::size_t nthOrder, int direction, class selector_tag>
class NthPartialDeriv; 

namespace internal{

// Evaluator  
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct Evaluator<fdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>> : public EvaluatorBase< fdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag> >
{
  const fdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>& m_xpr; 
  Evaluator(const fdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>& xpr): m_xpr(xpr){}
  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const fdm::Coordinate<N>& coords) const 
  {
    return Eigen::Map<const Eigen::Matrix<fdm::Scalar, Eigen::Dynamic, 1>>(weights + weights_per_order * _nthOrder, weights_per_order); 
  }
}; 


// Traits
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct traits_impl<fdm::linops::NthPartialDeriv<_nthOrder,_direction, selector_tag>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = 0; 
  static constexpr bool is_timedep = false; 
  static constexpr int direction = _direction; 
  static constexpr std::size_t maxOrder = _nthOrder; 
  typedef selector_tag node_selector_tag; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm

namespace Eigen{
namespace internal{

template<std::size_t _nthOrder, int _direction, class selector_tag>
struct traits<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>
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

template<std::size_t _nthOrder, int _direction, class selector_tag>
struct evaluator<fdm::linops::NthPartialDeriv<_nthOrder, _direction, selector_tag>> 
  : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>, 
  public evaluator_base<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>> 
{
  struct InnerIterator
    : public evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::InnerIterator
  {
    enum { CoeffReadCost = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::CoeffReadCost, Flags = evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::Flags };

    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::InnerIterator(eval, row_idx)
    {}
  }; 
  evaluator(const fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>& xpr_d)
    : evaluator<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>(xpr_d)
  {}
}; 

} // end namespace internal 
} // end namespace Eigen 

#include "PartialDerivBase.hpp"

namespace fdm{ 
namespace linops{

template<std::size_t _nthOrder, int _direction, class selector_tag = fdm::linops::Centered<0>>
class NthPartialDeriv : public PartialDerivBase<NthPartialDeriv<_nthOrder,_direction,selector_tag>>
{
  public:
    // Friends ----------------- 
    friend Eigen::internal::evaluator<NthPartialDeriv>; 
    friend fdm::linops::internal::EvaluatorBase<NthPartialDeriv>; 
    friend fdm::linops::internal::Evaluator<NthPartialDeriv>;

    // Member Data -------------- 
    static constexpr int direction = _direction; 
    static constexpr std::size_t order = _nthOrder; 
}; 

} // end namespace linops 
} // end namespace fdm 

#endif // NthPartialDeriv.hpp 