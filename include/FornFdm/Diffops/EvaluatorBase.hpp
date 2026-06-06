// Evaluator.hpp
//
// Composition class to efficiently write node expressions into CSRMatrix
//
// JAF 4/21/2026 

#ifndef FORNFDM_DIFFOPS_EVALUATORBASE_H
#define FORNFDM_DIFFOPS_EVALUATORBASE_H

#include<array>
#include "../utilities/Fornberg2.hpp"
#include "../utilities/FornbergStackCalc.hpp"
#include "../Types.hpp"
#include "../Traits.hpp"
#include "../Mesh.hpp" 
#include "../Coordinate.hpp"
#include "NodeSelector.hpp"

namespace fornfdm{
namespace linops{
namespace internal{

// declare as empty struct. specialized by individual types. should always inherit from EvaluatorBase<Xpr>  
template<class Xpr>
struct Evaluator{}; 

template<class Xpr>
struct EvaluatorBase
{
  using traits_t = fornfdm::linops::internal::traits<Xpr>; 
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
    typename fornfdm::utils::FornbergStackCalc<NodeSelector<typename traits_t::node_selector_tag, numNodesMin>::numNodesMax, traits_t::maxOrder> m_calc; 
    fornfdm::Coordinate<traits_t::max_num_args_called> m_coords; 
    fornfdm::Real m_time; 

    public:
    // constructor
    Row(const Evaluator<Xpr>& eval, const fornfdm::Mesh* m, std::size_t row_idx, fornfdm::Real t)
      : m_eval(eval), 
      m_nodes(m->getAxis(traits_t::direction), (row_idx / eval.m_xpr.m_prod_before) % m->sizeOfDim(traits_t::direction)), 
      m_calc(m_nodes.x_bar, m_nodes.nodeValues.cbegin(), std::next(m_nodes.nodeValues.cbegin(), m_nodes.numNodesUsed)), 
      m_coords(m, row_idx), 
      m_time(t)
    {}

    // Member Functions ======================================

    inline const std::size_t& size() const { return m_nodes.numNodesUsed; } 

    // returns indices the nodes were selected from. needs to be transformed by user of Row struct   
    const auto& columnIndices() const { return m_nodes.nodeIndices; } 
    
    // uses m_eval to map fornberg weights + coords into eigen expression  
    auto values() const { return m_eval.evalWeightsCoordsTime(m_calc.getArray().data(), m_nodes.numNodesUsed, m_coords, m_time); }

    // value of non zeros / row index offset 
    inline const std::size_t& valuePtrOffset() const { return m_nodes.nonZerosOffset; }  

    // function to return an eigen map of innerIndexPtr / valuePtr into Eigen::VectorXd for SIMD writing.
    template<typename U>
    auto mapToEigen(U* ptr) const { return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>>(ptr, m_nodes.numNodesUsed); }
  }; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#endif // Evaluator.hpp