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
// #include<OutsideSteps/All.hpp> 
// #include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 
#include<TExprs/All.hpp> 
// #include<Solvers/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp> 

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 
  
  // domain mesh 
  auto domain = LinOps::make_mesh(0.0,10.0,21); 

  // times mesh
  auto times = LinOps::make_mesh(0.0,4.0,20); 

  // ICs 
  BumpFunc b{.L = 3.0, .R=7.0, .c=5.0,  .h=2.0}; 
  std::vector<Eigen::VectorXd> sols{ LinOps::make_Discretization(domain, b).values() }; 

  // LHS in time 
  auto Ut = TExprs::NthTimeDeriv<1>{}; 

  // RHS in space 
  auto expr = 0.5 * LinOps::NthDerivOp(2) - 0.5 * LinOps::NthDerivOp(1); 
  expr.set_mesh(domain); 

  // Solving with explicit steps manually... 
  auto exec = TExprs::make_Executor(Ut); 
  auto it = times->cbegin(); 
  exec.pushTimeRange(it,++it); 
  exec.pushSolution(sols[0]); 

  auto end = times->cend(); 
  for(; it!=end; ++it)
  {
    exec.pushTime(*it); 
    exec.calculate(*std::prev(it)); 
    Eigen::VectorXd next_sol = exec.getInvCoeff() * expr.GetMat() * exec.getCurrentSolution() + exec.getRhsExpression();

    next_sol[0] = next_sol[next_sol.size()-1] = 0.0;
     
    sols.push_back(next_sol); 
    exec.pushSolution(std::move(next_sol)); 
  }
  
  print_mat(sols,"Solutions"); 
};