// tests_diffusion.cpp
//
// test on Ut = Uxx with IC= sin(x) on [0,pi]
// compared to the exact solution. 
//
// JAF 5/17/2026 

// workaround to allow expressions of different rows/cols to be added 
// before setMesh() sets their rows/cols to be equal
#define eigen_assert(x)
#include<FornFdm/All.hpp>
#include<FornFdm/Utilities/PrintVec.hpp>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <FornFdm/EigenFdmPlugin.hpp>

#include<iostream>
#include<iomanip>
#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

using namespace fornfdm; 

// Mesh Suite ---------------------------------------- 
TEST(DiffusionSuite, SinInitialConditionAllSolvers){

  // IO manip
  std::cout << std::setprecision(3); 
  constexpr double pi = 3.14159265385; 
  // Domain + Time  
  fornfdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fornfdm::linspaced(30,0.0,pi), 1), 
    .times = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(300,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  utils::print_vec(args.initialConditions[0],"Initial"); 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  auto left = osteps::DirichletBC(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  solvers::ExplicitSolver exp_solver(Ut,Uxx,std::tie(bcs)); 
  solvers::ImplicitSolver imp_solver(Ut,Uxx,std::tie(bcs)); 
  solvers::CrankNicolsonSolver cn_solver(Ut,Uxx,std::tie(bcs)); 

  auto test_solver = [&](auto& t_solver)
  {
    auto sol = t_solver.calculate(args, solvers::LastSaver{}); 
    utils::print_vec(sol,"Solution"); 
    for(auto idx = 0; idx<args.mesh->sizeOfDim(0); ++idx)
    {
      double exact_val = std::sin(args.mesh->getAxis(0)[idx]) * std::exp(- (*args.times)[args.times->size()-1]); 
      ASSERT_NEAR(sol[idx], exact_val, 0.025);
    }
  }; 

  test_solver(exp_solver); 
  test_solver(imp_solver); 
  test_solver(cn_solver); 

};