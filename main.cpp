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

using std::cout, std::endl;

using namespace fdm; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 
  
  // domain mesh 
  auto domain = make_Mesh1D(0.0,10.0,21); 

  // times mesh
  auto times = make_Mesh1D(0.0,4.0,20); 

  // ICs 
  utils::BumpFunc b{.L = 3.0, .R=7.0, .c=5.0,  .h=2.0}; 
  std::vector<Eigen::VectorXd> sols{ make_Discretization(domain, b).values() }; 

  // LHS in time 
  auto Ut = texprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto expr = 0.5 * linops::NthDerivOp<2>{} - 0.5 * linops::NthDerivOp<1>{}; 
  expr.setMesh(domain); 

  // Solving with explicit steps manually... 
  auto exec = texprs::make_Executor(Ut); 
  auto it = times->cbegin(); 
  exec.pushTimeRange(it,++it); 
  exec.pushSolution(sols[0]); 

  auto end = times->cend(); 
  Eigen::VectorXd next_sol; 

  for(; it!=end; ++it)
  {
    exec.pushTime(*it); 
    exec.calculate(*std::prev(it)); 
    next_sol = exec.getInvCoeff() * expr.asMatrix() * exec.getCurrentSolution() + exec.getRhsExpression();

    next_sol[0] = next_sol[next_sol.size()-1] = 0.0;
     
    sols.push_back(next_sol); 
    exec.pushSolution(std::move(next_sol)); 
  }
  
  utils::print_mat(sols,"Solutions"); 
};
