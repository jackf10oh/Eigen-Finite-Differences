// diffusion.cpp
//
// JAF 8/31/2026

#include<iostream>
#include<memory>
#include<fornfdm/all.hpp>
#include<fornfdm/utilities/BumpFunc.hpp>

using namespace fornfdm; 
using namespace fornfdm::linops; 

int main()
{
  // uniform spacing on [-1,1], 1D mesh
  auto m = std::make_shared<Mesh>(linspaced(30,0.0,1.0), 1);

  // Uxx + 0.5 Ux
  auto rhs = 0.1 * NthPartialDeriv<2,0>{} + NthPartialDeriv<1,0>{}; 
  rhs.setMesh(m);

  // uniform through time. 
  auto times = std::make_shared<const Vector>(m->getAxis(0));

  // Ut order 1 time derivative
  texprs::NthTimeDeriv<1> lhs; 

  // Initial condition is sin(pi*x)
  auto f = [](double x){ return std::sin(x * 3.14159); };
  auto ics = discretize(m,f);

  osteps::BCPair bcs(osteps::Dirichlet(0.0),osteps::Dirichlet(0.0));

  // Implicit Euler  
  solvers::ImplicitSolver solver(lhs,rhs,std::tie(bcs)); 
  solvers::SolverArgs args{m,times,{ics}};

  // prints
  solver.calculate(args,solvers::PrintSaver{});
};