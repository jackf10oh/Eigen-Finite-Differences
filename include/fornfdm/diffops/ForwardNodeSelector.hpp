// ForwardNodeSelector
// 
// selects nodes for forward scheme for 
// as many rows as possible 
// 
// JAF 5/11/2026 

#ifndef FORNFDM_DIFFOPS_FORWARDNODESELECTOR_H
#define FORNFDM_DIFFOPS_FORWARDNODESELECTOR_H 

#include<array>
#include "../types.hpp"
#include "NodeSelector.hpp"

namespace fornfdm{
namespace linops{

// tag type to take a guranteed minimum. outside of internal 
template<std::size_t numNodesMin=0>struct Forward{};

namespace internal{

template<std::size_t N, std::size_t M>
struct promote_node_selector_tags<Forward<N>,Forward<M>>
{
  constexpr static bool is_match = true;
  using type = Forward<std::max(N,M)>;
};

// forward declaration
template< std::size_t numNodesMin >
struct ForwardNodeSelector; 

// Specialize the NodeSelector base on tag
template<std::size_t numNodesMin01, std::size_t numNodesMin02>
struct NodeSelector<Forward<numNodesMin02>, numNodesMin01> : public ForwardNodeSelector<std::max(numNodesMin01,numNodesMin02)>
{
  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)  
    : ForwardNodeSelector<std::max(numNodesMin01,numNodesMin02)>(c, idx) 
  {}
}; 

// default node selector use centered stencil everywhere except at boundaries
template< std::size_t numNodesMin >
struct ForwardNodeSelector
{
  // Member Data ----------------------
  static constexpr std::size_t numNodesMax = numNodesMin; 
  static constexpr std::size_t numNodesUsed = numNodesMin; 
  std::array<Eigen::Index, numNodesMax> nodeIndices;
  Eigen::Map<const fornfdm::Vector> nodeValues;
  fornfdm::Scalar x_bar; 
  std::size_t nonZerosOffset; 

  // Constructor ============================
  template<class Container>
  ForwardNodeSelector(const Container& c, std::size_t idx)
    : nonZerosOffset(numNodesUsed * idx), nodeValues(nullptr,0)
  {
    // x_bar is just the ith entry into container c 
    x_bar = c[idx]; 

    // select the indices + nodes + offset
    auto s = c.size() - numNodesUsed;  
    if(idx < s){
      // top rows use a forward stencil 
      new (&nodeValues) Eigen::Map<const fornfdm::Vector>(c.data() + idx, numNodesUsed); 
      for(auto j=0; j<numNodesUsed; ++j){
        auto node_idx = idx + j;
        nodeIndices[j] = node_idx; 
      }
    }
    else {
      // bottom rows pack tightly against right end point
      new (&nodeValues) Eigen::Map<const fornfdm::Vector>(c.data() + s, numNodesUsed); 
      for(auto j=0; j<numNodesUsed; ++j){
        auto node_idx = s + j;
        nodeIndices[j] = node_idx; 
      }
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