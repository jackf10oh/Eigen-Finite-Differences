// internals.cpp
//
// JAF 8/31/2026

#include<iostream>
#include<memory>
#include<fornfdm/all.hpp>

using namespace fornfdm::linops; 
using std::endl, std::cout; 

int main()
{
  using fornfdm::Coordinate; 
  using fornfdm::linops::internal::Evaluator;
  double t{}; Coordinate<1> x{}; // dummies
  
  double d_01[] = {1, 0, -1, 1}; // weights from Fornberg
  double d_02[] = {0, 1, 0, -.5, 0, .5, -1, 2, -1};

  // 1.)
  NthPartialDeriv<1,0> Ux{};               
  auto lam_01 = internal::Evaluator<NthPartialDeriv<1,0>>(Ux,t).createReader(x);
  // [](const double* d, int n, int stride){ return d[(order 1)*stride + n]; }
  cout << lam_01(d_01, 0, 2) << ", " << lam_01(d_01, 1, 2) << "\n";
  // -1, 1

  // 2.)
  NthPartialDeriv<2,0> Uxx{};           
  auto lam_02 = internal::Evaluator<NthPartialDeriv<2,0>>(Uxx,t).createReader(x);
  // [](const double* d, int n, int stride){ return d[(order 2)*stride + n]; }
  cout << lam_02(d_02, 0, 3) << ", " 
       << lam_02(d_02, 1, 3) << ", " 
       << lam_02(d_02, 2, 3) << "\n";
  // -1, 2, -1

  // 3.)
  auto xpr = Uxx + Ux; 
  internal::Evaluator<std::decay_t<decltype(xpr)>> eval(xpr,t);
  auto lam_03 = eval.createReader(x);
  // [lam_01, lam_02](...){ return binop(lam_01(...),lam_02(...)); }
  cout << lam_03(d_02, 0, 3) << ", " 
       << lam_03(d_02, 1, 3) << ", " 
       << lam_03(d_02, 2, 3) << "\n";
  // -1/2+(-1) , 0+2 , 1/2+(-1)
};