// Evaluator.hpp
//
// Composition class to efficiently write node expressions into CSRMatrix
//
// JAF 4/21/2026 

#ifndef DIFFOP_EVALUATORBASE_H
#define DIFFOP_EVALUATORBASE_H 

#include "../Traits.hpp"

namespace fdm{
  namespace linops{

// Holds (x,y,z) coords in different dimmensions 
template< std::size_t max_dims >
struct Coordinate
{
  // Member Data ------------------------------
  std::array<fdm::Scalar, max_dims> values; 

  // Constructor ----------------------------
  Coordinate(const Mesh* m, std::size_t row_idx)
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

  // Member Functions ------------------------
  template<class Callable>
  fdm::Scalar apply(const Callable& c) const {
    constexpr std::size_t N = linops::traits::callable_traits<Callable>::num_args; 
    return apply_impl(c, std::make_index_sequence<N>{}); 
  }

  template<class Callable, std::size_t... idxs>
  fdm::Scalar apply_impl(const Callable& c, std::index_sequence<idxs...>) const { return c(values[idxs]...); }
};

namespace internal{

// declare as empty struct. specialized by individual types. should always inherit from EvaluatorBase<Xpr>  
template<class Xpr>
struct Evaluator{}; 

// default node selector use centered stencil everywhere except at boundaries
template< std::size_t numNodesMin >
struct CoreNodeSelector
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
  CoreNodeSelector(const Container& c, std::size_t idx)
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

// holds (x1, x2, ..., xn) in the same dimension 
template<class Xpr, std::size_t numNodesMin>
struct NodeSelector{}; 

// Specialize the NodeSelector base on tags 
struct centered_selector_tag {};
template<std::size_t numNodesMin>
struct NodeSelector<centered_selector_tag, numNodesMin> : public CoreNodeSelector<numNodesMin>
{
  // Constructor ============================
  template<class Container>
  NodeSelector(const Container& c, std::size_t idx)  
    : CoreNodeSelector<numNodesMin>(c, idx) 
  {}
}; 
// TODO? // struct forward_selector_tag {}; // struct backward_selector_tag {};

template<class Xpr>
struct EvaluatorBase
{
  using traits_t = fdm::linops::internal::traits<Xpr>; 
  static constexpr std::size_t numNodesMin = traits_t::maxOrder+1;

  // estimate is for compressed matrix. need to handle block diagonals outside of Evaluator. 
  template<class Container>
  static std::size_t nonZerosEstimate(const Container& c){ return NodeSelector<typename traits_t::node_selector_tag, numNodesMin>::sumNodesPerRow(c); } 
  
  // holds const Evaluator&, nodes, coords, and fornberg weights 
  class Row
  {
    private:
    // member data  
    const Evaluator<Xpr>& m_eval; 
    NodeSelector<typename traits_t::node_selector_tag, numNodesMin> m_nodes; 
    typename fdm::utils::FornArrayCalc<NodeSelector<typename traits_t::node_selector_tag, numNodesMin>::numNodesMax, traits_t::maxOrder> m_calc; 
    Coordinate<traits_t::max_num_args_called> m_coords; 

    public:
    // constructor
    Row(const Evaluator<Xpr>& eval, const Mesh* m, std::size_t row_idx)
      : m_eval(eval), 
      m_nodes(m->getAxis(traits_t::direction), (row_idx / m->sizesMiddleProduct(0,traits_t::direction)) % m->sizeOfDim(traits_t::direction)), 
      m_calc(m_nodes.x_bar, m_nodes.nodeValues.cbegin(), std::next(m_nodes.nodeValues.cbegin(), m_nodes.numNodesUsed)), 
      m_coords(m, row_idx)
    {}

    // Member Functions ======================================; 

    inline const std::size_t& size() const { return m_nodes.numNodesUsed; } 

    // returns indices the nodes were selected from. needs to be transformed by user of Row struct   
    const auto& columnIndices() const { return m_nodes.nodeIndices; } 
    
    // uses m_eval to map fornberg weights + coords into eigen expression  
    auto values() const { return m_eval.evaluateWeightsAndCoords(m_calc.getArray().data(), m_nodes.numNodesUsed, m_coords); }

    // value of non zeros / row index offset 
    inline const std::size_t& valuePtrOffset() const { return m_nodes.nonZerosOffset; }  

    // function to return an eigen map of innerIndexPtr / valuePtr into Eigen::VectorXd for SIMD writing.
    template<typename U>
    auto mapToEigen(U* ptr) const { return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>>(ptr, m_nodes.numNodesUsed); }
  }; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm 

#endif // Evaluator.hpp