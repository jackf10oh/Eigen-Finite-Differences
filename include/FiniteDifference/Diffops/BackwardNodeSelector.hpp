// BackwardNodeSelector.hpp
// 
// selects nodes for backward scheme 
// for as many rows as possible 
// 
// JAF 5/11/2026 

#ifndef FORNFDM_DIFFOPS_BACKWARDNODESELECTOR_H
#define FORNFDM_DIFFOPS_BACKWARDNODESELECTOR_H

namespace fornfdm{
namespace linops{

// tag type to take a guranteed minimum. outside of internal 
template<std::size_t numNodesMin=0>struct Backward{};

namespace internal{

// forward declaration
template< std::size_t numNodesMin >
struct BackwardNodeSelector; 

// Specialize the NodeSelector base on tag
template<std::size_t numNodesMin01, std::size_t numNodesMin02>
struct NodeSelector<Backward<numNodesMin02>, numNodesMin01> : public BackwardNodeSelector<std::max(numNodesMin01,numNodesMin02)>
{
  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)  
    : BackwardNodeSelector<std::max(numNodesMin01,numNodesMin02)>(c, idx) 
  {}
}; 

// default node selector use centered stencil everywhere except at boundaries
template< std::size_t numNodesMin >
struct BackwardNodeSelector
{
  // Member Data ----------------------
  static constexpr std::size_t numNodesMax = numNodesMin; 
  static constexpr std::size_t numNodesUsed = numNodesMin; 
  std::array<Eigen::Index, numNodesMax> nodeIndices;
  std::array<fornfdm::Scalar, numNodesMax> nodeValues;
  fornfdm::Scalar x_bar; 
  std::size_t nonZerosOffset; 

  // Constructor ============================
  template<class Container>
  BackwardNodeSelector(const Container& c, std::size_t idx)
    : nonZerosOffset(numNodesUsed * idx)
  {
    // x_bar is just the ith entry into container c 
    x_bar = c[idx]; 

    // select the indices + nodes + offset
    auto s = c.size() - numNodesUsed;  
    if(idx < numNodesUsed){
      // top rows pack as tightly as possible against left endpoint
      for(auto j=0; j<numNodesUsed; ++j){
        nodeIndices[j] = j; 
        nodeValues[j] = c[j];    
      }
    }
    else{
      // bottom rows use a backward stencil 
      auto j = numNodesUsed;
      do
      {
        auto node_idx = idx - numNodesMin + j; 
        --j; 
        nodeIndices[j] = node_idx; 
        nodeValues[j] = c[node_idx]; 
      } while (j!=0);
    }
  }

  // Member Functions ----------------------------------- 

  // sum of # of nodes in each row. 
  template<class Container>
  static std::size_t sumNodesPerRow(const Container& c)
  {
    // our case we have an exact formula 
    return c.size() * numNodesMax; 
    // in general we always need to have an exact formula for the case when 
    // there's a kronecker product that needs nnz * number of repeats 
  }
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm

#endif 