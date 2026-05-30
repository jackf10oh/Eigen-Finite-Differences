// tests_solvers.cpp 
// 
// tests interpolate() free function 
// in 1D, 2D, and 3D linear solutions 
// 
// JAF 5/18/2026 

// workaround to allow expressions of different rows/cols to be added 
// before setMesh() sets their rows/cols to be equal
#define eigen_assert(x)
#include<FornFdm/All.hpp>
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

// Solver Suite ------------------------------------------------- 
TEST(SolverSuite, OneDimExplicitSolver){
  constexpr double pi = 3.14159265385; 
  // Domain + Time  
  fornfdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fornfdm::linspaced(30,0.0,pi), 1), 
    .times = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(300,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  auto left = osteps::DirichletBC(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  solvers::ExplicitSolver exp_solver(Ut,Uxx,std::tie(bcs)); 
  
  auto sol = exp_solver.calculate(args, solvers::LastSaver{}); 
  for(auto idx = 0; idx<args.mesh->sizeOfDim(0); ++idx)
  {
    double exact_val = std::sin(args.mesh->getAxis(0)[idx]) * std::exp(- (*args.times)[args.times->size()-1]); 
    ASSERT_NEAR(sol[idx], exact_val, 0.025);
  }
};

TEST(SolverSuite, OneDimImplicitSolver){
  constexpr double pi = 3.14159265385; 
  // Domain + Time  
  fornfdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fornfdm::linspaced(30,0.0,pi), 1), 
    .times = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(300,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  auto left = osteps::DirichletBC(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  solvers::ImplicitSolver imp_solver(Ut,Uxx,std::tie(bcs)); 
  
  auto sol = imp_solver.calculate(args, solvers::LastSaver{}); 
  for(auto idx = 0; idx<args.mesh->sizeOfDim(0); ++idx)
  {
    double exact_val = std::sin(args.mesh->getAxis(0)[idx]) * std::exp(- (*args.times)[args.times->size()-1]); 
    ASSERT_NEAR(sol[idx], exact_val, 0.025);
  }
};

TEST(SolverSuite, OneDimCrankNicolsonSolver){
  constexpr double pi = 3.14159265385; 
  // Domain + Time  
  fornfdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fornfdm::linspaced(30,0.0,pi), 1), 
    .times = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(300,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  auto left = osteps::DirichletBC(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  solvers::CrankNicolsonSolver cn_solver(Ut,Uxx,std::tie(bcs)); 
  
  auto sol = cn_solver.calculate(args, solvers::LastSaver{}); 
  for(auto idx = 0; idx<args.mesh->sizeOfDim(0); ++idx)
  {
    double exact_val = std::sin(args.mesh->getAxis(0)[idx]) * std::exp(- (*args.times)[args.times->size()-1]); 
    ASSERT_NEAR(sol[idx], exact_val, 0.025);
  }
};

// Interpolation Suite ---------------------------------------- 
TEST(InterpolationSuite, HighDimLinearInterpolation){

  // 1D test 
  auto mesh_1d = fornfdm::make_Mesh(fornfdm::linspaced(21,-10.0,10.0),1); 
  auto mesh_1d_fine = fornfdm::make_Mesh(fornfdm::linspaced(200,-10.0,10.0),1); 
  auto lam_1d = [](double x){ return 10.54 + 32.7 * x; }; 
  auto vec_1d = make_Discretization(mesh_1d, lam_1d); 
  for(auto i = 0; i < mesh_1d_fine->sizeOfDim(0); ++i)
  {
    fornfdm::Coordinate<1> coord{mesh_1d_fine->getAxis(0)[i]};
    auto exact_val = coord.apply(lam_1d); 
    auto interp_val = fornfdm::solvers::interpolate(coord, vec_1d,mesh_1d.get());
    ASSERT_NEAR(interp_val, exact_val, 1e-4);
  }

  // 2D test 
  auto mesh_2d = fornfdm::make_Mesh(fornfdm::linspaced(21,-10.0,10.0),2); 
  auto mesh_2d_fine = fornfdm::make_Mesh(fornfdm::linspaced(200,-10.0,10.0),2); 
  auto lam_2d = [](double x, double y){ return 10.54 + 32.7 * x + 17.2 * y; }; 
  auto vec_2d = make_Discretization(mesh_2d, lam_2d); 
  for(auto i = 0; i < mesh_2d_fine->sizeOfDim(0); ++i)
  {
    for(auto j = 0; j < mesh_2d_fine->sizeOfDim(1); ++j)
    {
      fornfdm::Coordinate<2> coord{mesh_2d_fine->getAxis(0)[i], mesh_2d_fine->getAxis(1)[j]};
      auto exact_val = coord.apply(lam_2d); 
      auto interp_val = fornfdm::solvers::interpolate(coord, vec_2d,mesh_2d.get());
      ASSERT_NEAR(interp_val, exact_val, 1e-4);
    }
  }


  // 3D test
  auto mesh_3d = fornfdm::make_Mesh(fornfdm::linspaced(10,-10.0,10.0),3); 
  auto mesh_3d_fine = fornfdm::make_Mesh(fornfdm::linspaced(50,-10.0,10.0),3); 
  auto lam_3d = [](double x, double y, double z){ return 10.54 + 32.7 * x + 17.2 * y - 26.49 * z; }; 
  auto vec_3d = make_Discretization(mesh_3d, lam_3d); 
  for(auto i = 0; i < mesh_3d_fine->sizeOfDim(0); ++i)
  {
    for(auto j = 0; j < mesh_3d_fine->sizeOfDim(1); ++j)
    {
      for(auto k = 0; k < mesh_3d_fine->sizeOfDim(2); ++k)
      {
        fornfdm::Coordinate<3> coord{mesh_3d_fine->getAxis(0)[i], mesh_3d_fine->getAxis(1)[j], mesh_3d_fine->getAxis(2)[k]};
        auto exact_val = coord.apply(lam_3d); 
        auto interp_val = fornfdm::solvers::interpolate(coord, vec_3d,mesh_3d.get());
        ASSERT_NEAR(interp_val, exact_val, 1e-4);
      }
    }
  }

};

TEST(InterpolationSuite, OneDimInterpolator){

  constexpr double pi = 3.14159265385; 
  // Domain + Time  
  fornfdm::solvers::SolverArgs args{
    .mesh = std::make_shared<const fornfdm::Mesh>(fornfdm::linspaced(30,0.0,pi), 1), 
    .times = std::make_shared<const Eigen::VectorXd>(fornfdm::linspaced(300,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { std::move(v) }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  auto left = osteps::DirichletBC(0.0); 
  auto right = left;
  osteps::BCPair bcs(left,right); 

  solvers::CrankNicolsonSolver cn_solver(Ut,Uxx,std::tie(bcs)); 
  solvers::Interpolator interp(std::move(cn_solver), args);  
  const auto& axis = args.mesh->getAxis(0); 
  auto mesh_redone = fornfdm::make_Mesh(fornfdm::linspaced(25,axis[0], axis[axis.size()-1]), 1); 

  const auto& t_axis = *args.times; 
  auto times_redone = std::make_shared<const fornfdm::Vector>(fornfdm::linspaced(25,t_axis[0],t_axis[t_axis.size()-1])); 

  for(auto i = 0; i < mesh_redone->sizeOfDim(0); ++i)
  {
    for(auto j = 0; j < times_redone->size(); ++j)
    {
      double x = mesh_redone->getAxis(0)[i]; 
      double t = (*times_redone)[j]; 
      double exact_val = std::sin(x) * std::exp(-t); 
      double interpolated_val = interp.solAt<1>(t, {x}); 
      ASSERT_NEAR(interpolated_val, exact_val, 0.025);
    }
  }
};