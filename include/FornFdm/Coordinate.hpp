// Coordinate.hpp 
//
// Given an index from 0,1,2,...,N 
// where N = x_axis.size * ... * z_axis.size 
// packs coordinate of first N dimension 
// into an array 
//
// JAF 5/6/2026 

#ifndef FORNFDM_COORDINATE_H
#define FORNFDM_COORDINATE_H 

namespace fornfdm{

class Mesh; 

// Holds (x,y,z) coords in different dimmensions 
template< std::size_t numDimsMax >
struct Coordinate
{
  // Member Data ------------------------------
  std::array<fornfdm::Scalar, numDimsMax> values; 

  // Constructors ----------------------------

  // from Mesh + row_idx in flattened idx space  
  Coordinate(const Mesh* m, std::size_t row_idx);

  // from x,y,z values 
  template<typename... Xs, typename = std::enable_if_t<sizeof...(Xs)<=numDimsMax && (std::is_convertible_v<Xs, fornfdm::Scalar> && ...)> >
  Coordinate(Xs&&... xs)
    : values{ {std::forward<Xs>(xs)...} }
  {/* all args must be convertible to fornfdm::Scalar. sizeof...(Args) must be <= numDimsMax*/}

  // Member Functions ------------------------
  template<class Callable>
  fornfdm::Scalar apply(const Callable& c) const; 

  template<class Callable, class ArgType>
  fornfdm::Scalar applyBindFirst(const Callable& c, ArgType t) const; 

  private:
  // Implementations ------------------ 
  template<class Callable, std::size_t... idxs>
  fornfdm::Scalar apply_impl(const Callable& c, std::index_sequence<idxs...>) const; 

  template<class Callable, class ArgType, std::size_t... idxs>
  fornfdm::Scalar applyBindFirst_impl(const Callable& c, ArgType t, std::index_sequence<idxs...>) const; 
};

} // end namespace fornfdm 

#include "Mesh.hpp"

namespace fornfdm{ 

template<std::size_t numDimsMax>
Coordinate<numDimsMax>::Coordinate(const Mesh* m, std::size_t row_idx)
{
  if constexpr(numDimsMax > 0){
    std::size_t rolling_product = m->sizeOfDim(0); 
    values[0] = m->getAxis(0)[row_idx % rolling_product];  

    for(std::size_t ith_dim=1; ith_dim < numDimsMax; ++ith_dim){
      std::size_t s = m->sizeOfDim(ith_dim); 
      values[ith_dim] = m->getAxis(ith_dim)[(row_idx/rolling_product) % s];  
      rolling_product *= s; 
    }
  }
}

template<std::size_t numDimsMax>
template<class Callable>
fornfdm::Scalar Coordinate<numDimsMax>::apply(const Callable& c) const 
{
  constexpr std::size_t N = fornfdm::internal::callable_traits<Callable>::arity;
  if constexpr(N == 0){
    return c();
  }  
  else{
    return apply_impl(c, std::make_index_sequence<N>{}); 
  }
}

template<std::size_t numDimsMax>
template<class Callable, std::size_t... idxs>
fornfdm::Scalar Coordinate<numDimsMax>::apply_impl(const Callable& c, std::index_sequence<idxs...>) const 
{
  return c(values[idxs]...); 
}

template<std::size_t numDimsMax>
template<class Callable, class ArgType>
fornfdm::Scalar Coordinate<numDimsMax>::applyBindFirst(const Callable& c, ArgType t) const
{
  constexpr std::size_t N = fornfdm::internal::callable_traits<Callable>::arity;
  static_assert(N>0, "Must call applyBindFirst on callable F with arity >= 1"); 
  if constexpr(N==1){
    return c(t); 
  }
  else{
    return applyBindFirst_impl(c, t, std::make_index_sequence<N-1>{});
  }
}

template<std::size_t numDimsMax>
template<class Callable, class ArgType, std::size_t... idxs>
fornfdm::Scalar Coordinate<numDimsMax>::applyBindFirst_impl(const Callable& c, ArgType t, std::index_sequence<idxs...>) const
{
  return c(t,values[idxs]...); 
}

} // end namespace fornfdm 

#endif 