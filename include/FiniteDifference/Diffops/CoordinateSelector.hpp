// DiffopsUtils.hpp
//
// a couple of template declarations that are specialized by the Diffops library 
// TODO there's no use specializing this class based off of zero args in callables? 
// we can pass a default construct std::array to evaluateWeights with very little cost 
//
// JAF 4/18/2026 

#ifndef COORDINATESELECTOR_H
#define COORDINATESELECTOR_H

namespace fdm{
  namespace linops{
    namespace internal{ 
      
template<class U>
struct CoordinateSelector
{
  // Member Data 
  static constexpr std::size_t max_num_args_called = fdm::internal::traits<U>::max_num_args_called; 
  std::array<fdm::Scalar, max_num_args_called> coordinates; 

  // Constructor 
  CoordinateSelector(const Mesh* m, std::size_t row_idx)
  {

    std::size_t rolling_product = m->sizeOfDim(0); 
    coordinates[0] = m->getAxis(0)[row_idx % rolling_product];  

    for(std::size_t ith_dim=0; ith_dim<max_num_args_called; ++ith_dim){
      std::size_t s = m->sizeOfDim(ith_dim); 
      coordinates[ith_dim] = m->getAxis(ith_dim)[(row_idx/rolling_product) % s];  
      rolling_product *= s; 
    }
  }
};

    } // end namespace internal 
  } // end namespace linops 
} // end namespace fdm 

#endif // CoordinateSelector.hpp 