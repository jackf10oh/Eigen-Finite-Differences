// Fornberg.hpp
//
// header file for fornberg algorithm
//
// JAF 10/29/2025

#ifndef FORNBERGALGO_H
#define FORNBERGALGO_H

#include<cmath>
#include<vector>
#include<string> 
#include<iostream> 

// ------------------------------------------------------------------------
template<template<typename T,typename TAlloc> class Container, typename Real, typename RealAlloc>
std::vector<Real> LagrangeWeights(Real x_bar, const Container<Real,RealAlloc> nodes)
{
  // number of nodes 
  int N = nodes.size();
  if(N==0) throw std::runtime_error("must supply >= 1 nodes");

  // Vector of weights for each node (n=0, ..., N-1)
  auto Weights = std::vector<Real>(N, 0.0);

  // Using 1 nodes is just a flat line
  Weights[0] = 1; 

  // c4 will hold an accumulation of (x_bar - node[0]) * ... * (x_bar - node[n-1])
  Real c4=1; 
  // c5 will hold an accumulation of (node[n]-node[0]) * ... * (node[n]-node[n-1])
  Real c5; 

  // for number of nodes n=2, ..., N (first node was zero index)
  for(int n=1; n<N; n++)
  {
    // reset c5. it depends on node[n]
    c5=1; 

    // loop over all previous nodes 0, ..., n-1
    for(int old_n=0; old_n<n; old_n++)
    { 
      // all previous weights F(x_bar) are updated according to the rule 
      // F[n,v](x_bar) = ((x_bar-node[n]) / (node[v]-node[n])) * F[n-1,v](x)
      // (3.6) -------------------------------
      Weights[old_n] = ((x_bar-nodes[n])/(nodes[old_n]-nodes[n]))*Weights[old_n];
      
      // accumulation is updated as well 
      c5 *= (nodes[n]-nodes[old_n]);

      // next old node
    }

    // adding last node to accumulation
    c4 *= (x_bar-nodes[n-1]) ;
    // the newest weight at weight n is calculated
    // (3.2) ----------------------------
    Weights[n] = (c4/*product of all (x_bar - node)*/) / (c5/*product of all (nodes[old_n]-nodes[n])*/);

    // next new node
  }

  // Weights now contains LaGrange Interpolant Polynomials for 
  // nodes a0, a1, ..., an evaluated at x_bar
  return Weights; 
}

// ----------------------------------------------------------------------------------
template<template<typename T,typename TAlloc> class Container, typename Real, typename RealAlloc>
std::vector<Real> FornWeights(Real x_bar, int order, const Container<Real,RealAlloc> nodes) 
{
  if(order<0) throw std::runtime_error("derivative order must be >= 0");

  // number of nodes 
  int N = nodes.size();
  if(N==0) throw std::runtime_error("must supply >= 1 nodes");

  // Matrix of Order+1 rows, N cols
  // row m from Weights is the coeffs for derivative of order m (m= 0, ... , order)
  auto Weights = std::vector(order+1, std::vector<Real>(N, 0.0));

  // Using 1 nodes is just a flat line
  Weights[0][0] = 1; 

  // c1 holds old c2 for next loop
  Real c1=1.0; 
  // c2 will hold an accumulation of (node[n]-node[0]) * ... * (node[n]-node[n-1])
  Real c2; 
  // c3 holds the difference (nodes[new]-nodes[old])
  Real c3; 

  // for number of nodes n=2, ..., N (first node was zero index)
  for(int n=1; n<N; n++)
  {
    // std::cout << "N=" << n << std::endl;
    // c1 *= (x_bar-nodes[n-1]);
    // reset c5. it depends on node[n]
    c2=1.0; 
    // loop over all previous nodes 0, ..., n-1
    for(int old_n=0; old_n<n; old_n++)
    { 
      // accumulation is updated 
      c3 = (nodes[n]-nodes[old_n]); 
      c2 *= c3;

      // if old_n == n-1 we must use the very last old column 
      // to update the newest nodes column
      // this is because of the formulas. we cannot overwrite them yet!
      if(old_n==n-1)
      {
        // std::cout << "SPECIAL" << std::endl;
        // the newest weight at weight n is calculated
        // its derivates are calculated now too 
        for(int m=std::min(n,order); m>=1; m--)
        {
          // (3.9) ------------------------------------- 
          Weights[m][n] = (c1/c2) * (m*Weights[m-1][n-1] - (nodes[n-1]-x_bar)*Weights[m][n-1]); 
          // std::cout <<"old_n:"<<old_n+1<<" m:"<<m<<" weight[m][n]:" << Weights[m][n]<<std::endl;
        }
        // (3.7) ------------------------------
        Weights[0][n] = (c1/c2) * (x_bar - nodes[n-1]) * Weights[0][n-1];
        // std::cout <<"old_n:"<<old_n+1<<" m:0 weight[0][n]:" << Weights[0][n]<<std::endl;
        // std::cout << "END SPECIAL" << std::endl;      
      }

      // all previous weights F(x_bar) are updated according to the rule 
      // F[n,v](x_bar) = ((x_bar-node[n]) / (node[v]-node[n])) * F[n-1,v](x)
      // the derivates are updated as well 
      for(int m=std::min(n,order); m>=1; m--)
      {
        // (3.8) ----------------------------------------
        Weights[m][old_n] = (((nodes[n]-x_bar)*Weights[m][old_n]) - (m*Weights[m-1][old_n])) / c3;
        // std::cout <<"old_n:"<<old_n<<" m:"<<m<<" weight[m][old_n]:" << Weights[m][old_n]<<std::endl;
      }
      // (3.6) -------------------------------------------------
      Weights[0][old_n] = ((x_bar-nodes[n])/(nodes[old_n]-nodes[n]))*Weights[0][old_n];
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
  return Weights[order]; 
}

#endif 