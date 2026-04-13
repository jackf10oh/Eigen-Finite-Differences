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
#include<FiniteDifference/All.hpp> 

#include<FiniteDifference/Utilities/PrintVec.hpp> 
#include<FiniteDifference/Utilities/BumpFunc.hpp> 

#include<FiniteDifference/Solvers/CrankNicolsonSolver.hpp> 
#include<FiniteDifference/Solvers/Interpolator.hpp> 

using std::cout, std::endl;

using namespace fdm; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  fdm::solvers::SolverArgs args{
    .mesh = make_Mesh1D(0.0,10.0,40), 
    .times = make_Mesh1D(0.0, 4.0, 8)
  }; 

  // Initial Conditions  
  utils::BumpFunc b{.L = 4.0, .R=6.0, .c=5.0,  .h=1.0}; 
  auto v = make_Discretization(args.mesh, b).values(); 
  args.initialConditions = { v }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto expr = 0.2 * linops::NthDerivOp<2>{} - 0.5 * linops::NthDerivOp<1>{}; 

  // Boundary Conditions 
  auto bcs = osteps::BCPair(osteps::DirichletBC(0.0),osteps::DirichletBC(0.0)); 

  // Solving ...
  // solvers::ExplicitSolver my_solver(Ut,expr,std::tie(bcs)); 
  // solvers::ImplicitSolver my_solver(Ut,expr,std::tie(bcs)); 
  solvers::CrankNicolsonSolver my_solver(Ut,expr,std::tie(bcs)); 

  utils::print_vec(args.initialConditions[0],"ICs"); 
  auto sol = my_solver.calculate(args, solvers::LastSaver{}); 
  utils::print_vec(sol, "Sol"); 

  // my_solver.calculate(args, solvers::PrintSaver{}); 
  
  // auto time_taken = my_solver.calculate(args, solvers::TimerSaver{}); 
  // cout << "milliseconds: " << time_taken.count() << endl;  

  // std::size_t N = 40; 
  // double sum = 0; 
  // for(auto i=0; i<N; ++i) sum += my_solver.calculate(args, solvers::TimerSaver{}).count(); 
  // cout << "Average time: " << (sum/N) << " ms" << endl; 

  solvers::Interpolator my_interp(std::move(my_solver), args); 
  std::cout << "interp: ["; 
  auto it = args.mesh->cbegin(); 
  auto end = std::prev(args.mesh->cend()); 
  for(; it!=end;++it) cout << my_interp.SolAt(4.0, *it) << ", "; 
  cout << my_interp.SolAt(4.0, *it) << "]" << endl; ;

};
