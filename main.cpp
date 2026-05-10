// main.cpp
//
// JAF 12/8/2025

// workaround to allow expressions of different rows/cols to be added 
// before setMesh() sets their rows/cols to be equal
#define eigen_assert(x)
#include<FiniteDifference/All.hpp>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <FiniteDifference/EigenFdmPlugin.hpp> 

#include<iostream>
#include<iomanip>
#include<FiniteDifference/Utilities/PrintVec.hpp> 
#include<FiniteDifference/Utilities/BumpFunc.hpp>
#include<Eigen/SparseCore> // macro plugin takes effect. 

using namespace fdm; 

using std::endl, std::cout; 

int main()
{
 // iomanip 
  std::cout << std::setprecision(3); 
  
  // Domain + Time  
  fdm::solvers::SolverArgs args{
    .mesh = make_Mesh(fdm::linspaced(41,-5.0,5.0), 2), 
    .times = std::make_shared<const fdm::Vector>(fdm::linspaced(101,0.0,3.0))
  }; 

  // Initial Conditions  
  auto v = make_Discretization(args.mesh, 0.0); 
  args.initialConditions = { v, v }; 

  // LHS in time 
  auto Utt = texprs::NthTimeDeriv<2>{}; 

  // RHS in space 
  auto Uxx = linops::NthPartialDeriv<2,0>{}; 
  auto Uyy = linops::NthPartialDeriv<2,1>{}; 
  auto expr = Uxx + Uyy; 

  // Boundary Conditions 
  auto left = osteps::RobinBC(1.0,-1.0,0.0); 
  auto right = osteps::RobinBC(1.0,1.0,0.0);
  osteps::BCPair bc_pair(left,right); 
  osteps::BCList bcs(bc_pair,bc_pair); 

  // Forcing Terms 
  fdm::utils::BumpFunc bump{.L = -1.0, .R = 1.0, .c =0.0, .h = std::sqrt(5), .focus=10}; 
  osteps::ForcingTerm forcing = [bump](double t, double x, double y)
  {
    return std::sin(6.28318*t) * bump(x) * bump(y);
  }; 

  // Solving ...
  // solvers::ExplicitSolver my_solver(Utt,expr,std::tie(bcs)); 
  // solvers::ImplicitSolver my_solver(Utt,expr,std::tie(forcing, bcs)); 
  solvers::CrankNicolsonSolver my_solver(Utt,expr,std::tie(forcing, bcs)); 

  // 1D through time 
  // my_solver.calculate(args, solvers::PrintSaver{}); 

  // 1D Print 
  // auto sol = my_solver.calculate(args, solvers::LastSaver{}); 
  // utils::print_vec(args.initialConditions[0],"ICs"); 
  // utils::print_vec(sol, "Sol"); 

  // 2D print 
  auto sol = my_solver.calculate(args, solvers::LastSaver{}); 
  utils::print_mat(args.mesh->makeOneDimViews(args.initialConditions[0], 0), "Init"); 
  utils::print_mat(args.mesh->makeOneDimViews(sol, 0), "solution"); 
 
  // Time to last sol
  // auto time_taken = my_solver.calculate(args, solvers::TimerSaver{}); 
  // cout << "milliseconds: " << time_taken.count() << endl;  

  // Average time to last sol
  // std::size_t N = 40; 
  // double sum = 0; 
  // for(auto i=0; i<N; ++i) sum += my_solver.calculate(args, solvers::TimerSaver{}).count(); 
  // cout << "Average time: " << (sum/N) << " ms" << endl; 
};
