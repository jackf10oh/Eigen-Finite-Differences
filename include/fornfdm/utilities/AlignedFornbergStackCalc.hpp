// fornberg.hpp
//
// Redo of FornbergCalc that uses compile time fixed size array instead of vector
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 4/3/2026

#ifndef FORNFDM_UTILS_ALIGNEDFORNBERGSTACKCALC_H
#define FORNFDM_UTILS_ALIGNEDFORNBERGSTACKCALC_H

#include<array>
#include<cstdint>
#include<cmath>
#include<cassert>
#include<Eigen/Core>
#include "../types.hpp"
#include "aligned_fornberg.hpp"

namespace fornfdm{
  namespace utils{

// stateful Fornberg weight calculator. owns a fixed size array  
template<std::size_t M, std::size_t N>
class AlignedFornbergStackCalc
{
  public:
    // Member Data ----------------------------------------------------------
    static constexpr std::size_t alignment = EIGEN_MAX_ALIGN_BYTES / sizeof(fornfdm::Scalar);
    static constexpr std::size_t order = N;                                           // maximum order of derivative stencil
    static constexpr std::size_t numNodesMax = alignment*((M+alignment-1)/alignment); // first multiple of S such that k*S >= M.         
  private:
    EIGEN_ALIGN_MAX std::array<fornfdm::Scalar, numNodesMax*(N+1)> m_arr;  // single allocation of memory rows*cols big 
    std::size_t m_nodes_used;                                              // stores how many nodes were actually used in the algorithm
  
  public:
    // Constructors + Destructor =========================================================
    AlignedFornbergStackCalc()
    { 
      static_assert(M >= N + 1, "FornbergArrayCalc requires NUM_NODES >= ORDER + 1"); 
    }
    
    template<typename Iter>
    AlignedFornbergStackCalc(fornfdm::Scalar x_bar, Iter start, Iter end)
    { 
      static_assert(M >= N + 1, "FornbergArrayCalc requires NUM_NODES >= ORDER + 1"); 
      calculate(x_bar,start,end); 
    }
    
    AlignedFornbergStackCalc(const AlignedFornbergStackCalc& other)=default; 
    
    // destructor 
    ~AlignedFornbergStackCalc()=default; 
    
    // Member Funcs ======================================================================================
    
    // Const getter to stored weights 
    const auto& getArray() const { return m_arr; }
    auto getNumNodesUsed() const { return m_nodes_used; }  

    // Updates m_arr to contain weights up to order n
    template<typename Iter>
    void calculate(fornfdm::Scalar x_bar, Iter start, Iter end)
    {
      // make sure distance(start,end) <= numNodesMax 
      auto d = std::distance(start,end); 
      assert((d <= numNodesMax) && "AlignedFornbergStackCalc error: distance(start,end) > numNodesMax");  
      m_nodes_used = d; 
      fornfdm::utils::aligned_fornberg(start,end,x_bar,order,numNodesMax,m_arr.begin()); 
    }
};

  } // end namespace utils 
} // end namespace fornfdm 

#endif // FornbergArrayCalc.hpp