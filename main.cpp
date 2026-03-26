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
// #include<TExprs/All.hpp> 
// #include<Solvers/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

#include "OutsideSteps/include/OutsideSteps/StepContexts.hpp"
// #include "OutsideSteps/include/OutsideSteps/OStepBase.hpp" 
// #include "OutsideSteps/include/OutsideSteps/BoundaryConds1D/BCPair.hpp"
// #include "OutsideSteps/include/OutsideSteps/BoundaryConds1D/DirichletBC.hpp"   


using std::cout, std::endl;

int main()
{
  // iomanip 
  // std::cout << std::setprecision(3); 
  auto m = LinOps::make_mesh(0.0, 10.0, 11); 
  int x, r; 

  auto t01 = OSteps::make_time(); 
  auto t02 = OSteps::make_time(1.0); 
  auto t03 = OSteps::make_time(1.0,2.0); 
  auto t04 = OSteps::make_time(1.0,2.0, m); // creates its own owning copy of m ... 
  
  auto ctx01 = OSteps::make_context(); 
  auto ctx02 = OSteps::make_context(m); // creates its own owning copy of m ...  
  auto ctx03 = OSteps::make_context(m, &x); 
  auto ctx04 = OSteps::make_context(m, &x, &r);  

  print_vec(*t04.container, "time container"); 
  print_vec(*ctx04.getMesh(), "mesh domain"); 
};