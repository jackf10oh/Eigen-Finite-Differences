// FornbergCalc.hpp
//
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 10/29/2025

#ifndef FORNFDM_UTILS_FORNBERGCALC_H
#define FORNFDM_UTILS_FORNBERGCALC_H

#include<cmath>
#include "../types.hpp"
#include "fornberg.hpp"

namespace fornfdm{
  namespace utils{

// stateful Fornberg weight calculator. only allocate memory at creation 
class FornbergCalc
{
  public:
    // Member Data ----------------------------------------------------------
    std::size_t m_order;                // maximum order of derivative stencil  
    std::size_t m_n_nodes;              // number of nodes to use in approximation
    fornfdm::Vector m_vec;          // single allocation of memory rows*cols big 
  public:
    // Constructors + Destructor =========================================================
    FornbergCalc()=delete;
    FornbergCalc(std::size_t max_nodes, std::size_t max_order_init)
      : m_n_nodes(max_nodes), m_order(max_order_init), m_vec(m_n_nodes*(m_order+1))
    {}

    template<typename Iter>
    FornbergCalc(fornfdm::Scalar x_bar, Iter start, Iter end, std::size_t max_order_init)
    {
      calculate(x_bar, start, end, m_order); // updates m_n_nodes, m_order, etc...
    }

    // copy
    FornbergCalc(const FornbergCalc& other)=default; 

    // destructor 
    ~FornbergCalc()=default; 
    
    // Member Funcs ======================================================================================
    const auto& getVector() const {return m_vec; }
    auto order() const {return m_order; }
    auto getNumNodesUsed() const {return m_n_nodes; }

    // Updates m_arr to contain weights up to order n
    template<typename Iter>
    void calculate(fornfdm::Scalar x_bar, Iter start, Iter end, std::size_t order)
    {
      // make sure distance(start,end) <= numNodesMax 
      auto d = std::distance(start,end); 
      m_n_nodes = d; 
      m_order = order;
      m_vec.resize(d * (order + 1));
      fornfdm::utils::fornberg(start,end,x_bar,order,m_vec.begin()); 
    }
};

  } // end namespace utils 
} // end namespace fornfdm 

#endif // FornbergCalc.hpp