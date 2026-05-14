// Coordinate.hpp 
//
// Given an index from 0,1,2,...,N 
// where N = x_axis.size * ... * z_axis.size 
// packs coordinate of first N dimension 
// into an array 
//
// JAF 5/6/2026 

#ifndef FDM_COORDINATE_H
#define FDM_COORDINATE_H 

namespace fdm{

class Mesh; 

// Holds (x,y,z) coords in different dimmensions 
template< std::size_t numDimsMax >
struct Coordinate
{
  // Member Data ------------------------------
  std::array<fdm::Scalar, numDimsMax> values; 

  // Constructors ----------------------------

  // from Mesh + row_idx in flattened idx space  
  Coordinate(const Mesh* m, std::size_t row_idx);

  // from x,y,z values 
  template<typename... Xs>
  Coordinate(Xs&&... xs)
    : values{ {std::forward<Xs>(xs)...} }
  {
    static_assert(sizeof...(Xs)<=numDimsMax, "Must construct from <= numDimsMax in Coordinate"); 
    static_assert((std::is_convertible_v<Xs, fdm::Scalar> && ...), "All args must be convertible to fdm::Scalar"); 
  }

  // Member Functions ------------------------
  template<class Callable>
  fdm::Scalar apply(const Callable& c) const; 

  template<class Callable, std::size_t... idxs>
  fdm::Scalar apply_impl(const Callable& c, std::index_sequence<idxs...>) const; 
};

} // end namespace fdm 

#include "Mesh.hpp"

namespace fdm{ 

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
fdm::Scalar Coordinate<numDimsMax>::apply(const Callable& c) const 
{
  constexpr std::size_t N = fdm::internal::callable_traits<Callable>::arity;
  if constexpr(N == 0){
    return c();
  }  
  else{
    return apply_impl(c, std::make_index_sequence<N>{}); 
  }
}

template<std::size_t numDimsMax>
template<class Callable, std::size_t... idxs>
fdm::Scalar Coordinate<numDimsMax>::apply_impl(const Callable& c, std::index_sequence<idxs...>) const 
{
  return c(values[idxs]...); 
}

} // end namespace fdm 

#endif 