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
#include<Solvers/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // Defining Meshes + ICs 
  Solvers::SolverArgs args{
    .domain_mesh_ptr = LinOps::make_mesh(0.0, 10.0, 31), // start, end, nsteps 
    .time_mesh_ptr = LinOps::make_mesh(0.0,4.0, 21), 
    .ICs = {} 
  }; 
  // bump on [3,5] with maximum at (4, 3)
  BumpFunc f{.L = 3.0, .R = 5.0, .c=4.0, .h=10};
  args.ICs = { make_Discretization(args.domain_mesh_ptr, f).values() }; 
  
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
  Solvers::ImplicitSolver my_solver(time_expr, space_expr, std::tie(bcs));
  my_solver.Calculate(args, Solvers::PrintWrite{} ); 

  // Interpolating ---------------------------------------------------------
  Solvers::Interpolator my_interp(time_expr, space_expr, std::tie(bcs), args); 
  my_interp.FillVals(); 
  print_mat(my_interp.StoredData(), "Solutions through time");
  // std::cout << "solution at (t,x) :" << interp.SolAt(t,x) << std::endl; 
};