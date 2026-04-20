// NodeSelector.hpp 
//
// responsible for deciding which nodes to use in Fornberg algorithm. 
// currently uses a forward stencil in top rows, centered stencil in middle,
// and backward stencil in bottom rows
//
// JAF 4/18/2026 

#ifndef DIFFOPSNODESELECTOR_H
#define DIFFOPSNODESELECTOR_H 

#include<iostream>

namespace fdm{
  namespace linops{
    namespace internal{

template< std::size_t numNodesMin >
struct NodeSelector
{
  // Member Data ----------------------
  static constexpr std::size_t numNodesMax = 2*((numNodesMin)/2)+1; 
  std::size_t numNodesUsed; 
  std::array<Eigen::Index, numNodesMax> nodeIndices;
  std::array<fdm::Scalar, numNodesMax> nodeValues;
  double x_bar; 
  std::size_t nonZerosOffset; 

  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)
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