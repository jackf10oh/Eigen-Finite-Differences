// examples.cpp
//
// JAF 12/8/2025

#include<iostream>
#include<iomanip>
#include<fornfdm/plugin.hpp>
#include<fornfdm/all.hpp>
#include<fornfdm/utilities/print.hpp> 
#include<fornfdm/utilities/BumpFunc.hpp>
#include<Eigen/SparseCore> // macro plugin takes effect. 

using namespace fornfdm; 

using std::endl, std::cout; 

int main()
{
  // IO manip
  std::cout << std::setprecision(3); 

  // Domain + Times  
  constexpr double pi = 3.14159265385; 
  fornfdm::solvers::SolverArgs args{
    /*.mesh=*/ make_Mesh(fornfdm::linspaced(20,0.0,pi), 1), 
    /*.times=*/ std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(100,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = fornfdm::discretize(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 
  Uxx.setMesh(args.mesh); 

  fornfdm::CSRMatrix stencil = 0.5 * Uxx.evalTime(20.0); 

  cout << stencil << endl;

  cout << Uxx << endl;

};
