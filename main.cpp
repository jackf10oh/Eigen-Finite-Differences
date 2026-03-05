// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include<LinOps/All.hpp> 
#include<OutsideSteps/All.hpp> 
#include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 
#include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

#include<Solvers/ExplicitSolver.hpp> 
#include<Solvers/ImplicitSolver.hpp> 
#include<Solvers/Interpolator.hpp>

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0;
  int n_gridpoints = 31;

  // mesh in space [0, 10] 
  auto my_mesh = LinOps::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time [0, 2]
  auto time_mesh = LinOps::make_mesh(0.0, 2.0, 17); 

  // Initializing IC discretizations -------------------------------------------------------
  // bump on [3,5] with maximum at (4, 3)
  BumpFunc f{.L = 3.0, .R = 5.0, .c=4.0, .h=10};  

  // solution U( t=0 ) 
  auto U0 = LinOps::make_Discretization(my_mesh, f).values();
  
  // Solver Args expect std::vec of ICs for first N time steps 
  std::vector<Eigen::VectorXd> ics = {std::move(U0)};  

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(1); 

  // building RHS expression -----------------------------------------------------
  using D = LinOps::NthDerivOp;
  auto space_expr = 0.2 * D(2) - 1.0 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  auto left = OSteps::DirichletBC(0.0); 
  auto right = OSteps::DirichletBC(0.0); 
  auto bcs = OSteps::BCPair(left,right); 

  // Solving --------------------------------------------------------------------- 
  SolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .ICs = ics, 
  }; 

  ImplicitSolver my_solver(time_expr, space_expr, std::tie(bcs));
  my_solver.Calculate(args, PrintWrite{} ); 

  // Interpolating ---------------------------------------------------------
  Interpolator my_interp(time_expr, space_expr, std::tie(bcs), args); 
  my_interp.FillVals(); 
  print_mat(my_interp.StoredData(), "Solutions through time");
  // std::cout << "solution at (t,x) :" << interp.SolAt(t,x) << std::endl; 
};