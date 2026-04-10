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

  auto m = linops::make_Mesh1D(0.0,10.0,11); 
  auto m2 = linops::make_MeshXD(0.0,4.0,5,2); 
  
  auto f1 = [](){ return 2.0; }; 
  auto f2 = [](double t){ return t; }; 
  auto f3 = [](double t, double x){ return t + x*x; }; 
  
  // linops::TimeDepCoeff c01 = f1; 
  linops::TimeDepCoeff c02 = f2; 
  linops::TimeDepCoeff c03 = f3; 

  linops::IOp I; 
  
  auto expr = c02 * I; 

  cout << "Identity is time dep? " << I.isTimeDep << endl; 
  cout << "TimeDepCoeff is time dep? " << c02.isTimeDep << endl;
  cout << "expr is time dep? " << expr.isTimeDep << endl;  
};