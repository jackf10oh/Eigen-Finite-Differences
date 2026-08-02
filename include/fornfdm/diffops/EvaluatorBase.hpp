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
  const fornfdm::Real m_time;
  // Constructors + Destructor ------- 
  EvaluatorBase()=delete;
  EvaluatorBase(fornfdm::Real t)
    : m_time(t)
  {}
  EvaluatorBase(const EvaluatorBase& other)=default;
  ~EvaluatorBase()=default;

  using traits_t = fornfdm::linops::internal::traits<Xpr>; 
  static constexpr std::size_t numNodesMin = traits_t::max_order+1;

  // estimate is for compressed matrix. need to handle block diagonals outside of Evaluator. 
  template<class Container>
  static std::size_t totalNonZeros(const Container& c){ return NodeSelector<typename traits_t::node_selector_tag, numNodesMin>::sumNodesPerRow(c); } 
  
  // holds const Evaluator&, nodes, coords, and fornberg weights 
  class Row
  {
    private:
    // types
    using Selector = NodeSelector<typename traits_t::node_selector_tag, numNodesMin>;  
    using FornCalc = typename fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::max_order>;
    // member data  
    const Evaluator<Xpr>& m_eval; 
    Selector m_nodes; 
    FornCalc m_calc; 
    fornfdm::Coordinate<traits_t::max_arity> m_coords; 
    decltype (m_eval.createReader(m_coords)) m_reader;

    public:
    // Constructors -------------------- 
    [[deprecated("use 2nd constructor. no longer performing  (idx / prod_before) % axis_size from idx. new constructor takes (eval, mesh*, node_idx, row_idx, time) where node_idx is idx into 1D axis and row_idx is index into matrix.")]]
    Row(const Evaluator<Xpr>& eval, const fornfdm::Mesh* m, std::size_t row_idx)
      : m_eval(eval), 
      m_nodes(m->getAxis(traits_t::direction), (row_idx / eval.m_xpr.m_prod_before) % m->sizeOfDim(traits_t::direction)),
      m_calc(m_nodes.x_bar, m_nodes.nodeValues.cbegin(), std::next(m_nodes.nodeValues.cbegin(), m_nodes.numNodesUsed)), 
      m_coords(m, row_idx), 
      m_reader(m_eval.createReader(m_coords))
    {}

    Row(const Evaluator<Xpr>& eval, const fornfdm::Mesh* m, std::size_t node_idx, std::size_t row_idx)
      : m_eval(eval), 
      m_nodes(m->getAxis(traits_t::direction), node_idx),
      m_calc(m_nodes.x_bar, m_nodes.nodeValues.cbegin(), std::next(m_nodes.nodeValues.cbegin(), m_nodes.numNodesUsed)), 
      m_coords(m, row_idx), 
      m_reader(m_eval.createReader(m_coords))
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