// NthPartialDeriv.hpp 
//
// Concrete class for Nth Partial derivative in direction. 
// 
// JAF 4/17/2026 

#ifndef FORNFDM_DIFFOPS_NTHPARTIALDERIV_H
#define FORNFDM_DIFFOPS_NTHPARTIALDERIV_H

#include<cstdint>
#include "../types.hpp"
#include "traits.hpp"
#include "EvaluatorBase.hpp"
#include "PartialDerivBase.hpp"
#include "EigenEvaluatorImpl.hpp"
#include "CenteredNodeSelector.hpp"

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
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  { 
    return [](const fornfdm::Scalar* weights, std::size_t idx, std::size_t stride){ return weights[_nthOrder*stride + idx]; };
  }

  struct ExactReader
  {    
    template<std::size_t... orders>
    fornfdm::Scalar operator()(const fornfdm::Scalar* data, std::size_t idx, std::size_t stride, std::index_sequence<orders...>) const
    {
      constexpr std::size_t offset = locate_order<_nthOrder, std::index_sequence<orders...>>::value; 
      return data[offset*stride + idx];
    }   
  };

  template<std::size_t N>
  ExactReader createExactReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return ExactReader{};
  }
}; 

// traits
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct traits_impl<fornfdm::linops::NthPartialDeriv<_nthOrder,_direction, selector_tag>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = 0; 
  static constexpr bool is_timedep = false; 
  static constexpr int direction = _direction; 
  static constexpr std::size_t max_order = _nthOrder;
  typedef std::index_sequence<_nthOrder> orders; 
  typedef selector_tag node_selector_tag; 
}; 

// map_to_base
template<std::size_t _nthOrder, int _direction, class selector_tag>
struct map_to_base_tag<
  fornfdm::linops::NthPartialDeriv<_nthOrder,_direction, selector_tag>, 
  std::enable_if_t<(_direction == 0)>
>{ using type = LeftKroneckerTag; }; 

template<std::size_t _nthOrder, int _direction, class selector_tag>
struct map_to_base_tag<
  fornfdm::linops::NthPartialDeriv<_nthOrder,_direction, selector_tag>, 
  std::enable_if_t<(_direction != 0)>
>{ using type = DoubleKroneckerTag; }; 

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
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>
{
  using XprType = fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>; 
  using Impl = typename fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::NthPartialDeriv<_nthOrder, _direction,selector_tag>>; 
  using InnerIterator = typename Impl::InnerIterator; 
  evaluator(const XprType& xpr)
    : Impl(xpr)
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
    // Member Data -------------- 
    static constexpr int direction = _direction; 
    static constexpr std::size_t order = _nthOrder; 
}; 

} // end namespace linops 
} // end namespace fornfdm 

#endif // NthPartialDeriv.hpp 