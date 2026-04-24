// NthPartialDeriv.hpp 
//
// Concrete class for Nth Partial derivative in direction. 
// 
// JAF 4/17/2026 

#ifndef NTHPARTIALDERIV_H
#define NTHPARTIALDERIV_H

#include "EvaluatorBase.hpp"
#include "../Traits.hpp"

namespace fdm{
namespace linops{

  // Forward Declaration ------------------------------------------------- 
template<std::size_t nthOrder, int direction>
class NthPartialDeriv; 

namespace internal{

// Evaluator  
template<std::size_t _nthOrder, int _direction>
struct Evaluator<fdm::linops::NthPartialDeriv<_nthOrder,_direction>> : public EvaluatorBase< fdm::linops::NthPartialDeriv<_nthOrder,_direction> >
{
  Evaluator(const fdm::linops::NthPartialDeriv<_nthOrder,_direction>& xpr){}
  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const fdm::linops::Coordinate<N>& coords) const 
  {
    return Eigen::Map<const Eigen::Matrix<fdm::Scalar, Eigen::Dynamic, 1>>(weights + weights_per_order * _nthOrder, weights_per_order); 
  }
}; 


// Traits
template<std::size_t _nthOrder, int _direction>
struct traits_impl<fdm::linops::NthPartialDeriv<_nthOrder,_direction>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = 0; 
  static constexpr bool is_timedep = false; 
  static constexpr int direction = _direction; 
  static constexpr std::size_t maxOrder = _nthOrder; 
  typedef centered_selector_tag node_selector_tag; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm

namespace Eigen{
namespace internal{

template<std::size_t _nthOrder, int _direction>
struct traits<fdm::linops::NthPartialDeriv<_nthOrder, _direction>>
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
    Flags = Eigen::RowMajor | NestByRefBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

} // end namespace internal 
} // end namespace Eigen 

#include "PartialDerivBase.hpp"

namespace fdm{ 
namespace linops{

template<std::size_t nthOrder, int direction>
class NthPartialDeriv : public PartialDerivBase<NthPartialDeriv<nthOrder,direction>>
{}; 

} // end namespace linops 
} // end namespace fdm 

#endif // NthPartialDeriv.hpp 