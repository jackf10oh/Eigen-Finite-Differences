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

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

#include<LinOps/All.hpp> 

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  auto m = linops::make_mesh(0.0,10.0,11); 
  auto m2 = linops::make_meshXD(0.0,4.0,5,2); 
  
  auto f1 = [](){ return 1.0; }; 
  auto f2 = [](double x){ return x; }; 
  auto f3 = [](double x){ return x*x; }; 
  
  linops::AutonomousCoeff c01 = f1; 
  linops::AutonomousCoeff c02 = f2; 
  linops::AutonomousCoeff c03 = f3; 

  // c01.setMesh1D(m); 
  // cout << c01.asMatrix() << endl; 

};