// NthPartialDeriv.hpp 
//
// Concrete class for Nth Partial derivative in direction. 
// 
// JAF 4/17/2026 

#ifndef FORNFDM_DIFFOPS_NTHPARTIALDERIV_H
#define FORNFDM_DIFFOPS_NTHPARTIALDERIV_H

#include "EvaluatorBase.hpp"
#include "CenteredNodeSelector.hpp"
#include "../Traits.hpp"

namespace fornfdm{
namespace linops{

  // Forward Declaration ------------------------------------------------- 
template<std::size_t nthOrder, int direction, class selector_tag>
class NthPartialDeriv; 

namespace internal{

// Evaluator  
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct Evaluator<fornfdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>> : public EvaluatorBase< fornfdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag> >
{
  const fornfdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>& m_xpr; 
  Evaluator(const fornfdm::linops::NthPartialDeriv<_nthOrder,_direction,selector_tag>& xpr): m_xpr(xpr){}
  template<std::size_t N>
  auto evalWeightsCoordsTime(const fornfdm::Scalar* weights, std::size_t weights_per_order, const fornfdm::Coordinate<N>& coords, fornfdm::Real t) const 
  {
    return Eigen::Map<const Eigen::Matrix<fornfdm::Scalar, Eigen::Dynamic, 1>>(weights + weights_per_order * _nthOrder, weights_per_order); 
  }
}; 


// Traits
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct traits_impl<fornfdm::linops::NthPartialDeriv<_nthOrder,_direction, selector_tag>>
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
} // end namespace fornfdm

namespace Eigen{
namespace internal{

template<std::size_t _nthOrder, int _direction, class selector_tag>
struct traits<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>
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

template<std::size_t _nthOrder, int _direction, class selector_tag>
struct evaluator<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction, selector_tag>> 
  : public evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>, 
  public evaluator_base<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>> 
{
  struct InnerIterator
    : public evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::InnerIterator
  {
    enum { CoeffReadCost = evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::CoeffReadCost, Flags = evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::Flags };

    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>::InnerIterator(eval, row_idx)
    {}
  }; 
  evaluator(const fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>& xpr_d)
    : evaluator<fornfdm::linops::PartialDerivBase<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>>(xpr_d)
  {}
}; 

} // end namespace internal 
} // end namespace Eigen 

#include "PartialDerivBase.hpp"

namespace fornfdm{ 
namespace linops{

template<std::size_t _nthOrder, int _direction, class selector_tag = fornfdm::linops::Centered<0>>
class NthPartialDeriv : public PartialDerivBase<NthPartialDeriv<_nthOrder,_direction,selector_tag>>
{
  public:
    // Friends ----------------- 
    friend Eigen::internal::evaluator<NthPartialDeriv>; 
    friend fornfdm::linops::internal::EvaluatorBase<NthPartialDeriv>; 
    friend fornfdm::linops::internal::Evaluator<NthPartialDeriv>;

    // Member Data -------------- 
    static constexpr int direction = _direction; 
    static constexpr std::size_t order = _nthOrder; 
}; 

} // end namespace linops 
} // end namespace fornfdm 

#endif // NthPartialDeriv.hpp 