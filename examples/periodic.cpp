// periodic.cpp
//
// JAF 8/31/2026

#include<iostream>
#include<memory>
#include<fornfdm/all.hpp>
#include<fornfdm/utilities/BumpFunc.hpp>

using namespace fornfdm; 
using namespace fornfdm::linops; 

// Periodic with length = 2
struct Period{ static constexpr double value = 2.0; };

int main()
{
  // uniform spacing on [0,2], 1D mesh
  auto m = std::make_shared<Mesh>(linspaced(40,0.0,1.95), 1);

  // Uxx + 0.5 Ux
  using S = Periodic<Period,0>;
  auto rhs = 0.005 * NthPartialDeriv<2,0,S>{} + 0.8 * NthPartialDeriv<1,0,S>{}; 

  // uniform through time. 
  auto times = std::make_shared<const Vector>(linspaced(201,0.0,10.0));

  // Ut order 1 time derivative
  texprs::NthTimeDeriv<1> lhs{}; 

  // Initial condition is a smooth bump centered at x=1.0
  utils::BumpFunc f{.L=0.5, .R=1.5, .c=1.0, .h=4.0}; 
  auto ics = discretize(m,f);

  // Implicit Euler  
  solvers::ImplicitSolver solver(lhs,rhs,std::tie()); 
  solvers::SolverArgs args{m,times,{ics}};

  // prints
  solver.calculate(args,solvers::PrintSaver{});
};