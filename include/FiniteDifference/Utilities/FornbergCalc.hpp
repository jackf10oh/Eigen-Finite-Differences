// Fornberg.hpp
//
// header file for stateful fornberg calc that only allocates on construction 
//
// JAF 10/29/2025

#ifndef FORNBERGCALC_H
#define FORNBERGCALC_H

#include<cmath>
#include<vector>
#include<string> 
#include<iostream> 
#include<Eigen/Dense> // Eigen::Map 

// stateful Fornberg weight calculator. only allocate memory at creation 
class FornCalc
{
  public:
    // Member Data ----------------------------------------------------------
    std::size_t m_order;                // maximum order of derivative stencil  
    std::size_t m_n_nodes;              // number of nodes to use in approximation
    std::vector<double> m_arr;          // single allocation of memory rows*cols big 
  public:
    // Constructors + Destructor =========================================================
    FornCalc()=delete;
    FornCalc(std::size_t max_nodes, std::size_t max_order_init=1)
      : m_n_nodes(max_nodes), m_order(max_order_init), m_arr(m_n_nodes*(m_order+1))
    {
      // if(m_order+1>max_nodes) throw std::invalid_argument("# Nodes must be > Order deriv in Fornberg algo"); 
      // for(double& val : m_arr) val=0.0; 
    };
    FornCalc(const FornCalc& other)=delete; 
    // destructor 
    ~FornCalc()=default; 
    
    // Member Funcs ======================================================================================
    template<typename Input_Iter>
    Eigen::Map<const Eigen::VectorXd> GetWeights(double x_bar, Input_Iter start, Input_Iter end, std::size_t order=1)
    {
      // Update m_arr weight values 
      Calculate(x_bar, start, end, order); 

      // return non owning view of last row of Weights;
      return Eigen::Map<const Eigen::VectorXd>(m_arr.data()+(m_n_nodes*order), m_n_nodes); 
    };

    // Updates m_arr to contain weights up to order n
    template<typename Input_Iter>
    void Calculate(double x_bar, Input_Iter start, Input_Iter end, std::size_t order=1)
    {
      // Matrix of Order+1 rows, N cols
      // row m from Weights is the coeffs for derivative of order m (m= 0, ... , order)
      // ---> we need atleast Order+1 * N entries in m_arr 
      if(std::distance(start,end)*(order+1) > m_arr.size()) throw std::runtime_error("nodes*(order+1) exceeds size of stored array!"); 
      // m_arr.resize(m_n_nodes*(order+1));

      // utility lambdas convert (i,j) -> index in flattened m_arr 
      auto entryRef = [this](std::size_t i, std::size_t j)->double&{ return this->m_arr[i*this->m_n_nodes+j];}; 

      // number of nodes 
      m_n_nodes = std::distance(start,end); 
      m_order = order;
      auto nodeRef = [&start](std::size_t i)-> const double& {return *(start+i);};
      
      // zero all stored entries
      for(auto& entry : m_arr) entry=0.0; 

      // Using 1 nodes is just a flat line
      entryRef(0,0) = 1; 

      // c1 holds old c2 for next loop
      double c1=1.0; 
      // c2 will hold an accumulation of (node[n]-node[0]) * ... * (node[n]-node[n-1])
      double c2; 
      // c3 holds the difference (nodes[new]-nodes[old])
      double c3; 

      // for number of nodes n=2, ..., N (first node was zero index)
      for(std::size_t n=1; n<m_n_nodes; n++)
      {
        // std::cout << "N=" << n << std::endl;
        // c1 *= (x_bar-nodes[n-1]);
        // reset c5. it depends on node[n]
        c2=1.0; 
        // loop over all previous nodes 0, ..., n-1
        for(std::size_t old_n=0; old_n<n; old_n++)
        { 
          // accumulation is updated 
          c3 = (nodeRef(n)-nodeRef(old_n)); 
          c2 *= c3;

          // if old_n == n-1 we must use the very last old column 
          // to update the newest nodes column
          // this is because of the formulas. we cannot overwrite them yet!
          if(old_n==n-1)
          {
            // std::cout << "SPECIAL" << std::endl;
            // the newest weight at weight n is calculated
            // its derivates are calculated now too 
            for(std::size_t m=std::min(n,m_order); m>=1; m--)
            {
              // (3.9) ------------------------------------- 
              entryRef(m,n) = (c1/c2) * ( m*entryRef(m-1,n-1) - (nodeRef(n-1)-x_bar)*entryRef(m,n-1) ); 
              // std::cout <<"old_n:"<<old_n+1<<" m:"<<m<<" weight[m][n]:" << Weights[m][n]<<std::endl;
            }
            // (3.7) ------------------------------
            entryRef(0,n) = (c1/c2) * (x_bar - nodeRef(n-1)) * entryRef(0,n-1);
            // std::cout <<"old_n:"<<old_n+1<<" m:0 weight[0][n]:" << Weights[0][n]<<std::endl;
            // std::cout << "END SPECIAL" << std::endl;      
          }

        // all previous weights F(x_bar) are updated according to the rule 
        // F[n,v](x_bar) = ((x_bar-node[n]) / (node[v]-node[n])) * F[n-1,v](x)
        // the derivates are updated as well 
        for(std::size_t m=std::min(n,order); m>=1; m--)
        {
          // (3.8) ----------------------------------------
          entryRef(m,old_n) = (((nodeRef(n)-x_bar)*entryRef(m,old_n)) - (m*entryRef(m-1,old_n))) / c3;
          // std::cout <<"old_n:"<<old_n<<" m:"<<m<<" weight[m][old_n]:" << Weights[m][old_n]<<std::endl;
        }
        // (3.6) -------------------------------------------------
        entryRef(0,old_n) = ((x_bar-nodeRef(n))/(nodeRef(old_n)-nodeRef(n)))*entryRef(0,old_n);
        // std::cout <<"old_n:"<<old_n<<" m:0 weight[0][old_n]:" << Weights[0][old_n] <<std::endl;

        // next old node
        }

        // accumulator is stored for next loop 
        c1 = c2; 
      
        // next new node
        // std::cout << std::endl;
      }

      // Weights now contains LaGrange Interpolant Polynomials for 
      // nodes a0, a1, ..., an evaluated at x_bar
      
    } // end Calculate ------------
};

#endif // FornbergCalc.hpp