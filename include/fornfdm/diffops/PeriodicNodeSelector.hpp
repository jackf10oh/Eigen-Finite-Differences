// PeriodicNodeSelector.hpp
// 
// Uses the centered scheme for all rows. 
// top/bottom rows wrap around and treat 
// the values in c as if they are on [0,Period::value]
// 
// i.e. c = [0,1,2,3,...8,9] with Period::value=10
// the top row may choose [0,1,-1], indices [0,1,9]
// with bottom row [10,8,9], indices [10,8,9] 
// 
// JAF 5/11/2026 

#ifndef FORNFDM_DIFFOPS_PERIODICNODESELECTOR_H
#define FORNFDM_DIFFOPS_PERIODICNODESELECTOR_H

#include<array>
#include "../types.hpp"
#include "NodeSelector.hpp"

namespace fornfdm{
namespace linops{

// tag type
// takes class with ::value=L to define a period on [0.0, L] 
// along with a minimum # of nodes
template<class Period, std::size_t numNodesMin=0>struct Periodic
{
  static constexpr decltype(Period::value) period = Period::value;
};

namespace internal{

template<class P, std::size_t m, class Q, std::size_t n>
struct promote_node_selector_tags<Periodic<P,m>,Periodic<Q,n>>
{
  constexpr static bool is_match = (P::value == Q::value);
  using type = Periodic<P, std::max(m,n)>;
};

// forward declaration
template<class P, std::size_t n>
struct PeriodicNodeSelector; 

// Specialize the NodeSelector base on tag
template<class P, std::size_t m, std::size_t n>
struct NodeSelector<Periodic<P, m>, n> : public PeriodicNodeSelector<P,std::max(m,n)>
{
  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)  
    : PeriodicNodeSelector<P,std::max(m,n)>(c, idx) 
  {}
}; 

// selects centered scheme in middle rows. top/bottom rows wrap around
template<class Period, std::size_t numNodesMin>
struct PeriodicNodeSelector
{
  // Member Data ----------------------
  static constexpr std::size_t numNodesMax = 2*((numNodesMin)/2)+1; 
  static constexpr std::size_t numNodesUsed = numNodesMax; 
  std::array<std::size_t, numNodesMax> nodeIndices;
  std::array<fornfdm::Scalar, numNodesMax> nodeValues;
  fornfdm::Scalar x_bar; 
  std::size_t nonZerosOffset; 

  // Constructor ============================
  template<class Container>
  PeriodicNodeSelector(const Container& c, std::size_t idx)
    : nonZerosOffset(idx * numNodesUsed)
  {
    // x_bar is just the ith entry into container c 
    x_bar = c[idx];

    // select the indices + nodes + offset
    constexpr std::size_t centered_skirt = (numNodesMin)/2; 
    if(idx<centered_skirt){
      // wraps around to end of c 
      
      // auto k=c.size()-idx
      std::size_t m = 1 + centered_skirt + idx;
      std::size_t j=0;
      for(; j < m; ++j){
        nodeIndices[j] = j;
        nodeValues[j] = c[j];
      }
      for(auto k = m; k<numNodesUsed; ++k){
        std::size_t node_idx = c.size() - numNodesUsed + k; 
        nodeIndices[j] = node_idx;
        nodeValues[j] = c[node_idx] - Period::value;
        ++j;
      }
    }
    else if(idx < c.size()-centered_skirt){
      // middle rows use a centered stencil 
      for(auto j=0; j<numNodesUsed; ++j){
        auto node_idx = idx-centered_skirt+j; 
        nodeIndices[j] = node_idx;
        nodeValues[j] = c[node_idx]; 
      }
    }
    else{
      // wraps around to beginning of c


      // will always have centered_skirt + 1 nodes available. 
      // starting from last idx. wraps to idx 0,1,2,...
      // idx == size() - 1 -> centered_skirt wrap arounds. 
      // idx == size() - 2 -> centered_skirt - 1 wrap arounds.
      // idx == size() - 3 -> centered_skirt - 2 wrap arounds.
      // ...
      std::size_t wraps = centered_skirt - (c.size() - 1 - idx);
      std::size_t k = 0;
      for(; k<wraps; ++k)
      {
        nodeIndices[k] = k;
        nodeValues[k] = Period::value + c[k];
      }
      std::size_t node_idx = idx - centered_skirt;
      for(; k<numNodesUsed; ++k)
      {
        nodeIndices[k] = node_idx;
        nodeValues[k] = c[node_idx];
        ++node_idx;
      }
    }
  }

  // Member Functions ----------------------------------- 

  // sum of # of nodes in each row. 
  template<class Container>
  static std::size_t sumNodesPerRow(const Container& c)
  {
    // our case we have an exact formula 
    constexpr std::size_t centered_skirt = (numNodesMin)/2;  
    return c.size() * (1+2*centered_skirt); 
    // in general we always need to have an exact formula for the case when 
    // there's a kronecker product that needs nnz * number of repeats 
  }
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm

#endif 