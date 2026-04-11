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
#include<FiniteDifference/All.hpp> 

#include<FiniteDifference/OutsideSteps/All.hpp> 

using std::cout, std::endl;

using namespace fdm; 

int main()
{
  // iomanip 
  // std::cout << std::setprecision(3); 
  auto m = fdm::make_MeshXD(0.0, 3.0, 5, 2); 

  auto mesh = fdm::make_Mesh1D(0.0,10.0, 11);
  
  fdm::utils::print_vec(*mesh, "mesh"); 


};
