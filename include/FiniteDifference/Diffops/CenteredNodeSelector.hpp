// CenteredNodeSelector.hpp
//
// Selects centered stencil as many rows 
// possible in the middle of the contianer. 
//
// JAF 5/11/2026 

#ifndef FDM_DIFFOPS_CENTEREDNODESELECTOR_H
#define FDM_DIFFOPS_CENTEREDNODESELECTOR_H 

#include "NodeSelector.hpp"
#include "EvaluatorBase.hpp"

namespace fdm{
namespace linops{

// tag type to take a guranteed minimum. outside of internal 
template<std::size_t numNodesMin=0> struct Centered{};

namespace internal{

// forward declaration
template< std::size_t numNodesMin >
struct CenteredNodeSelector; 

// Specialize the NodeSelector base on tag
template<std::size_t numNodesMin01, std::size_t numNodesMin02>
struct NodeSelector<Centered<numNodesMin02>, numNodesMin01> : public CenteredNodeSelector<std::max(numNodesMin01, numNodesMin02)>
{
  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)  
    : CenteredNodeSelector<std::max(numNodesMin01, numNodesMin02)>(c, idx) 
  {}
}; 

// default node selector use centered stencil everywhere except at boundaries
template< std::size_t numNodesMin >
struct CenteredNodeSelector
{
  // Member Data ----------------------
  static constexpr std::size_t numNodesMax = 2*((numNodesMin)/2)+1; 
  std::size_t numNodesUsed; 
  std::array<Eigen::Index, numNodesMax> nodeIndices;
  std::array<fdm::Scalar, numNodesMax> nodeValues;
  fdm::Scalar x_bar; 
  std::size_t nonZerosOffset; 

  // Constructor ============================
  template<class Container>
  CenteredNodeSelector(const Container& c, std::size_t idx)
  {
    constexpr std::size_t centered_skirt = (numNodesMin)/2; 
    // calculate how many nodes were actually used in the array. 
    numNodesUsed = (centered_skirt <= idx && idx < c.size()-centered_skirt) ? (1 + 2 * centered_skirt) : (numNodesMin); 

    // x_bar is just the ith entry into container c 
    x_bar = c[idx]; 

    // select the indices + nodes + offset
    if(idx<centered_skirt){
      nonZerosOffset = idx*numNodesUsed; 
      // top rows use a forward stencil 
      for(auto j=0; j<numNodesUsed; ++j){
        auto node_idx = idx + j;
        nodeIndices[j] = node_idx; 
        nodeValues[j] = c[node_idx];    
      }
    }
    else if(idx < c.size()-centered_skirt){
      nonZerosOffset = centered_skirt*numNodesMin + (idx-centered_skirt)*(1+2*centered_skirt); 
      // middle rows use a centered stencil 
      for(auto j=0; j<numNodesUsed; ++j){
        auto node_idx = idx-centered_skirt+j; 
        nodeIndices[j] = node_idx; 
        nodeValues[j] = c[node_idx]; 
      }
    }
    else{
      nonZerosOffset = centered_skirt*numNodesMin + (c.size()-2*centered_skirt)*(1+2*centered_skirt) + (idx-centered_skirt-(c.size()-2*centered_skirt))*(numNodesMin); 
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
    constexpr std::size_t centered_skirt = (numNodesMin)/2;  
    return 2*centered_skirt*(numNodesMin) + (c.size()-2*centered_skirt)*(1+2*centered_skirt); 
    // in general we always need to have an exact formula for the case when 
    // there's a kronecker product that needs nnz * number of repeats 
  }
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm

#endif 