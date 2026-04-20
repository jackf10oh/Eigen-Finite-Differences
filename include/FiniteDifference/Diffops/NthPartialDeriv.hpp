// NthPartialDeriv.hpp 
//
// Concrete class for Nth Partial derivative in direction. 
// 
// JAF 4/17/2026 

#ifndef NTHPARTIALDERIV_H
#define NTHPARTIALDERIV_H

namespace fdm{
namespace linops{

// Forward Declaration ------------------------------------------------- 
template<std::size_t nthOrder, int direction>
class NthPartialDeriv; 

} // end namespace linops 

namespace internal{

// Traits ------------------------------------------------------------
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
}; 

} // end namespace internal 
} // end namespace fdm

namespace Eigen{
namespace internal{

template<std::size_t nthOrder, int direction>
struct traits<fdm::linops::NthPartialDeriv<nthOrder,direction>> : traits<fdm::linops::PartialDerivBase<fdm::linops::NthPartialDeriv<nthOrder,direction>>>{}; 

} // end namespace internal 
} // end namespace Eigen 

#include "PartialDerivBase.hpp"

namespace fdm{ 
namespace linops{

template<std::size_t nthOrder, int direction>
class NthPartialDeriv : public PartialDerivBase<NthPartialDeriv<nthOrder,direction>>
{
  friend PartialDerivBase<NthPartialDeriv<nthOrder,direction>>; 
  // TODO NthPartialDeriv has a template parameter for NodeSelector and declares it inside the class.... 
  public:
    template<std::size_t numNodesMax, std::size_t numCoordsMax=0>
    auto evaluateWeightsAndCoords(
      const std::array<double, numNodesMax>& weights, 
      std::size_t weights_per_order, 
      const std::array<double,numCoordsMax>& coords={}) const 
    {
      using Mapped = Eigen::Map<const Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>>;  
      return Mapped(weights.data() + nthOrder * weights_per_order, weights_per_order); 
    }
}; 

} // end namespace linops 
} // end namespace fdm 

#endif // NthPartialDeriv.hpp 