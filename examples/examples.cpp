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
  fornfdm::solvers::SolverArgs args{
    /*.mesh=*/ make_Mesh(fornfdm::linspaced(20,0.0,3.14159), 1), 
    /*.times=*/ std::make_shared<const fornfdm::RealVector>(fornfdm::linspaced(100,0.0,0.5))
  }; 

  // Initial Conditions  
  auto v = fornfdm::discretize(args.mesh, [](double x){ return std::sin(x); }); 
  args.initialConditions = { v }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0,fornfdm::linops::Centered<5>>{}; 

  // Boundary Conditions 
  osteps::Dirichlet left(0.0), right(0.0); 
  osteps::BCPair bcs(left,right); 

  // Solving ...
  // solvers::ExplicitSolver my_solver(Ut,Uxx,std::tie(bcs)); 
  // solvers::ImplicitSolver my_solver(Utt,expr,std::tie(forcing, bcs)); 
  solvers::CrankNicolsonSolver my_solver(Ut,Uxx,std::tie(bcs)); 

  // 1D through time 
  my_solver.calculate(args, solvers::PrintSaver{}); 

  // 1D Print 
  // auto sol = my_solver.calculate(args, solvers::LastSaver{}); 
  // utils::print_vec(args.initialConditions[0],"Initial"); 
  // utils::print_vec(sol, "Sol"); 

  // 2D print 
  // auto sol = my_solver.calculate(args, solvers::LastSaver{}); 
  // utils::print_mat(args.mesh->makeOneDimViews(args.initialConditions[0], 0), "Init"); 
  // utils::print_mat(args.mesh->makeOneDimViews(sol, 0), "solution"); 
 
  // Time to last sol
  // auto time_taken = my_solver.calculate(args, solvers::TimerSaver{}); 
  // cout << "milliseconds: " << time_taken.count() << endl;  

  // Average time to last sol
  // std::size_t N = 40; 
  // double sum = 0; 
  // for(auto i=0; i<N; ++i) sum += my_solver.calculate(args, solvers::TimerSaver{}).count(); 
  // cout << "Average time: " << (sum/N) << " ms" << endl; 
};
