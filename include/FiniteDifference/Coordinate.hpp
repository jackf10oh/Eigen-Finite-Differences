// Coordinate.hpp 
//
//
//
// JAF 5/6/2026 

#ifndef FDM_COORDINATE_H
#define FDM_COORDINATE_H 

namespace fdm{

class Mesh; 

// Holds (x,y,z) coords in different dimmensions 
template< std::size_t max_dims >
struct Coordinate
{
  // Member Data ------------------------------
  std::array<fdm::Scalar, max_dims> values; 

  // Constructor ----------------------------
  Coordinate(const Mesh* m, std::size_t row_idx);

  // Member Functions ------------------------
  template<class Callable>
  fdm::Scalar apply(const Callable& c) const; 

  template<class Callable, std::size_t... idxs>
  fdm::Scalar apply_impl(const Callable& c, std::index_sequence<idxs...>) const; 
};

} // end namespace fdm 

#include "Mesh.hpp"

namespace fdm{ 

template<std::size_t max_dims>
Coordinate<max_dims>::Coordinate(const Mesh* m, std::size_t row_idx)
{
  if constexpr(max_dims > 0){
    std::size_t rolling_product = m->sizeOfDim(0); 
    values[0] = m->getAxis(0)[row_idx % rolling_product];  

    for(std::size_t ith_dim=1; ith_dim < max_dims; ++ith_dim){
      std::size_t s = m->sizeOfDim(ith_dim); 
      values[ith_dim] = m->getAxis(ith_dim)[(row_idx/rolling_product) % s];  
      rolling_product *= s; 
    }
  }
}

template<std::size_t max_dims>
template<class Callable>
fdm::Scalar Coordinate<max_dims>::apply(const Callable& c) const 
{
  constexpr std::size_t N = fdm::internal::callable_traits<Callable>::num_args;
  if constexpr(N == 0){
    return c();
  }  
  else{
    return apply_impl(c, std::make_index_sequence<N>{}); 
  }
}

template<std::size_t max_dims>
template<class Callable, std::size_t... idxs>
fdm::Scalar Coordinate<max_dims>::apply_impl(const Callable& c, std::index_sequence<idxs...>) const 
{
  return c(values[idxs]...); 
}



} // end namespace fdm 

#endif 