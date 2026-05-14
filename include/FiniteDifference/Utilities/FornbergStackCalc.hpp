// Fornberg.hpp
//
// Redo of FornbergCalc that uses compile time fixed size array instead of vector
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 4/3/2026

#ifndef FDM_UTILS_FORNBERGSTACKCALC_H
#define FDM_UTILS_FORNBERGSTACKCALC_H

#include<cmath>
#include<cassert>
#include<array>
#include<string> 
#include<iostream> 
#include "Fornberg2.hpp" // actual implementation of algorithm 

namespace fdm{
  namespace utils{

// stateful Fornberg weight calculator. owns a fixed size array  
template<std::size_t M, std::size_t N>
class FornbergStackCalc
{
  public:
    // Member Data ----------------------------------------------------------
    static constexpr std::size_t order = N;                // maximum order of derivative stencil  
    static constexpr std::size_t numNodesMax = M;              // number of nodes to use in approximation
  private:
    std::array<fdm::Scalar, M*(N+1)> m_arr;          // single allocation of memory rows*cols big 
    std::size_t m_nodes_used; // stores how many nodes were actually used in the algorithm
  public:
    // Constructors + Destructor =========================================================
    FornbergStackCalc()
    { 
      static_assert(M >= N + 1, "FornbergArrayCalc requires NUM_NODES >= ORDER + 1"); 
    }
    
    template<typename Iter>
    FornbergStackCalc(fdm::Scalar x_bar, Iter start, Iter end)
    { 
      static_assert(M >= N + 1, "FornbergArrayCalc requires NUM_NODES >= ORDER + 1"); 
      calculate(x_bar,start,end); 
    }
    
    FornbergStackCalc(const FornbergStackCalc& other)=default; 
    
    // destructor 
    ~FornbergStackCalc()=default; 
    
    // Member Funcs ======================================================================================
    
    // Const getter to stored weights 
    const auto& getArray() const { return m_arr; }
    auto getNumNodesUsed() const { return m_nodes_used; }  

    // Updates m_arr to contain weights up to order n
    template<typename Iter>
    void calculate(fdm::Scalar x_bar, Iter start, Iter end)
    {
      // make sure distance(start,end) <= numNodesMax 
      auto d = std::distance(start,end); 
      assert((d <= numNodesMax) && "FornbergStackCalc error: distance(start,end) > numNodesMax");  
      m_nodes_used = d; 
      fdm::utils::fornberg2(start,end,x_bar,order,m_arr.begin()); 
    }
};

  } // end namespace utils 
} // end namespace fdm 

#endif // FornbergArrayCalc.hpp