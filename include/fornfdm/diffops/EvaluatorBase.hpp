// Evaluator.hpp
//
// Composition class to efficiently write node expressions into CSRMatrix
//
// JAF 4/21/2026 

#ifndef FORNFDM_DIFFOPS_EVALUATORBASE_H
#define FORNFDM_DIFFOPS_EVALUATORBASE_H

#include<cstdint>
#include<array>
#include "../types.hpp"
#include "../traits.hpp"
#include "../Mesh.hpp" 
#include "../Coordinate.hpp"
#include "../utilities/fornberg.hpp"
#include "../utilities/FornbergStackCalc.hpp"
#include "traits.hpp"
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
  static std::size_t totalNonZeros(const Container& c){ return NodeSelector<typename traits_t::node_selector_tag, numNodesMin>::sumNodesPerRow(c); } 
  
  // holds const Evaluator&, nodes, coords, and fornberg weights 
  class Row
  {
    private:
    // types
    using Selector = NodeSelector<typename traits_t::node_selector_tag, numNodesMin>;  
    using FornCalc = typename fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::maxOrder>;
    // member data  
    const Evaluator<Xpr>& m_eval; 
    Selector m_nodes; 
    FornCalc m_calc; 
    fornfdm::Coordinate<traits_t::max_num_args_called> m_coords; 
    fornfdm::Real m_time; 
    decltype (m_eval.createReader(m_coords,m_time)) m_reader;

    public:
    // constructor
    Row(const Evaluator<Xpr>& eval, const fornfdm::Mesh* m, std::size_t row_idx, fornfdm::Real t)
      : m_eval(eval), 
      m_nodes(m->getAxis(traits_t::direction), (row_idx / eval.m_xpr.m_prod_before) % m->sizeOfDim(traits_t::direction)), 
      m_calc(m_nodes.x_bar, m_nodes.nodeValues.cbegin(), std::next(m_nodes.nodeValues.cbegin(), m_nodes.numNodesUsed)), 
      m_coords(m, row_idx), 
      m_time(t),
      m_reader(m_eval.createReader(m_coords,m_time))
    {}

    // Member Functions ======================================

    inline std::size_t size() const { return m_nodes.numNodesUsed; } 

    // returns indices the nodes were selected from. needs to be transformed by user of Row struct   
    inline std::size_t index(std::size_t ith) const { return m_nodes.nodeIndices[ith]; } 

    // uses m_reader to get value of a weight
    inline fornfdm::Scalar value(std::size_t ith) const { return m_reader(m_calc.getArray().data(), ith, m_nodes.numNodesUsed); }

    // value of non zeros / row index offset 
    inline std::size_t offset() const { return m_nodes.nonZerosOffset; }  
  }; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#endif // Evaluator.hpp