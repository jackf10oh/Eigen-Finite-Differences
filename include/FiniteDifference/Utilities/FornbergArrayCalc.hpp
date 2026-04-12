// Fornberg.hpp
//
// Redo of FornbergCalc that uses compile time fixed size array instead of vector
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 4/3/2026

#ifndef FORNBERGARRAYCALC_H
#define FORNBERGARRAYCALC_H

#include<cmath>
#include<array>
#include<string> 
#include<iostream> 
// #include<Eigen/Dense> // Eigen::Map 

namespace fdm{
  namespace utils{

// stateful Fornberg weight calculator. owns a fixed size array  
template<std::size_t M, std::size_t N>
class FornArrayCalc
{
  public:
    // Member Data ----------------------------------------------------------
    static constexpr std::size_t order = N;                // maximum order of derivative stencil  
    static constexpr std::size_t numNodes = M;              // number of nodes to use in approximation
  private:
    std::array<double, M*(N+1)> m_arr;          // single allocation of memory rows*cols big 
  public:
    // Constructors + Destructor =========================================================
    FornArrayCalc(){ static_assert(M+1 >= N, "FornbergArrayCalc requires NUM_NODES + 1 >= ORDER"); };
    FornArrayCalc(const FornArrayCalc& other)=default; 
    // destructor 
    ~FornArrayCalc()=default; 
    
    // Member Funcs ======================================================================================
    
    // Const getter to stored weights 
    const auto& getArray() const { return m_arr; };  

    // Updates m_arr to contain weights up to order n
    template<typename Input_Iter>
    void calculate(double x_bar, Input_Iter start, Input_Iter end)
    {
      // make sure distance(start,end) == numNodes 
      if(std::distance(start,end) != numNodes) throw std::runtime_error("FornbergArrayCalc error: distance(start,end) != numNodes");  

      // utility lambdas convert (i,j) -> index in flattened m_arr 
      auto entryRef = [this](std::size_t i, std::size_t j)->double& { return this->m_arr[i*numNodes+j];}; 
      auto nodeRef = [&start](std::size_t i)-> const double& {return *std::next(start,i); };
      
      // zero all stored entries
      for(auto& entry : m_arr) entry=0.0; 

      // Using 1 nodes is just a flat line
      entryRef(0,0) = 1; 

      // if using >= 2 nodes 
      if constexpr(numNodes >= 2)
      {
        // update first row according to legrend polynomials by hand  
        entryRef(0,0) = (x_bar - nodeRef(1)) / (nodeRef(0) - nodeRef(1)); 
        entryRef(0,1) = (x_bar - nodeRef(0)) / (nodeRef(1) - nodeRef(0));

        // derivative of legrange polynomial is just constant for 2 nodes i.e. straight line y = mx + b
        entryRef(1,0) = (1 / (nodeRef(0) - nodeRef(1)));
        entryRef(1,1) = (1 / (nodeRef(1) - nodeRef(0)));
      }

      // c1 holds old c2 for next loop
      double c1 = (nodeRef(1)-nodeRef(0))*(nodeRef(1)-nodeRef(0)); 
      // c2 will hold an accumulation of (node[n]-node[0]) * ... * (node[n]-node[n-1])
      double c2; 
      // c3 holds the difference (nodes[new]-nodes[old])
      double c3; 
    
      // for number of nodes n=3, ..., N (first node was zero index)
      if constexpr(numNodes >= 3){
        for(std::size_t n=2; n < numNodes; n++)
        {
          // c1 *= (x_bar-nodes[n-1]);
          // reset c5. it depends on node[n]
          c2=1.0; 
          // loop over all previous nodes 0, ..., n-1
          std::size_t old_n = 0; 
          // for(std::size_t old_n=0; old_n<n; old_n++)
          do
          { 
            // accumulation is updated 
            c3 = (nodeRef(n)-nodeRef(old_n)); 
            c2 *= c3;

            // all previous weights F(x_bar) are updated according to the rule 
            // F[n,v](x_bar) = ((x_bar-node[n]) / (node[v]-node[n])) * F[n-1,v](x)
            // the derivates are updated as well 
            for(std::size_t m=std::min(n,order); m>=1; m--)
            {
              // (3.8) ----------------------------------------
              entryRef(m,old_n) = (((nodeRef(n)-x_bar)*entryRef(m,old_n)) - (m*entryRef(m-1,old_n))) / c3;
            }
            // (3.6) -------------------------------------------------
            entryRef(0,old_n) = ((x_bar-nodeRef(n))/(nodeRef(old_n)-nodeRef(n)))*entryRef(0,old_n);

            // next old node
            ++old_n; 

          } while(old_n != n-1); 

          // accumulation is updated 
          c3 = (nodeRef(n)-nodeRef(old_n)); 
          c2 *= c3;

          // if old_n == n-1 we must use the very last old column 
          // to update the newest nodes column
          // this is because of the formulas. we cannot overwrite them yet!
          // the newest weight at weight n is calculated
          // its derivates are calculated now too 
          for(std::size_t m=std::min(n,order); m>=1; m--)
          {
            // (3.9) ------------------------------------- 
            entryRef(m,n) = (c1/c2) * ( m*entryRef(m-1,n-1) - (nodeRef(n-1)-x_bar)*entryRef(m,n-1) ); 
          }
          // (3.7) ------------------------------
          entryRef(0,n) = (c1/c2) * (x_bar - nodeRef(n-1)) * entryRef(0,n-1); 

          // finally update weights in 2nd to last column
          // the derivates are updated as well 
          for(std::size_t m=std::min(n,order); m>=1; m--)
          {
            // (3.8) ----------------------------------------
            entryRef(m,old_n) = (((nodeRef(n)-x_bar)*entryRef(m,old_n)) - (m*entryRef(m-1,old_n))) / c3;
          }
          // (3.6) -------------------------------------------------
          entryRef(0,old_n) = ((x_bar-nodeRef(n))/(nodeRef(old_n)-nodeRef(n)))*entryRef(0,old_n);

          // accumulator is stored for next loop 
          c1 = c2; 
        
          // next new node
        }
      }

      // Weights now contains LaGrange Interpolant Polynomials for 
      // nodes a0, a1, ..., an evaluated at x_bar
      
    } // end calculate ------------
};

  } // end namespace utils 
} // end namespace fdm 

#endif // FornbergArrayCalc.hpp