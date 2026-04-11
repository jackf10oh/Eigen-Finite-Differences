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

#include<FiniteDifference/Utilities/PrintVec.hpp> 
#include<FiniteDifference/Utilities/BumpFunc.hpp> 
#include<FiniteDifference/Utilities/RowMajorIdentityExpr.hpp> 
#include<FiniteDifference/Mesh1D.hpp> 

// #include<LinOps/All.hpp> 

using std::cout, std::endl;

using namespace fdm; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  auto mesh = make_Mesh1D(0.0,10.0, 11);
  
  print_vec(*mesh, "mesh"); 

  // auto tmp01 = linops::IOp(); 
  // auto tmp02 = 3.0 * tmp01; 
  // auto tmp03 = 3.0 * tmp02; 
  // auto rhs1 = RandOp().compose(tmp01); 
  // auto rhs2 = tmp03 - IOp().compose(rhs1); 

  // auto my_expr = rhs2; 

  // my_expr.setMesh(mesh); 
  // cout << my_expr.asMatrix() << endl; 
};