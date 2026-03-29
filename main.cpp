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
#include<OutsideSteps/All.hpp> 
// #include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 
// #include<TExprs/All.hpp> 
// #include<Solvers/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

int main()
{
  // iomanip 
  // std::cout << std::setprecision(3); 
  auto m = LinOps::make_meshXD(0.0, 3.0, 5, 2); 

  OSteps::DirichletBC bc01(5.0); 
  OSteps::DirichletBC bc02(7.0); 

  OSteps::BCPair p01(bc01, bc01); 
  OSteps::BCPair p02(bc02, bc02);

  OSteps::BCList bcs(p01,p02);  

  LinOps::MatrixStorage_t M = LinOps::DirectionalNthDerivOp(m,1,0).GetMat(); 
  bcs.MatBeforeStep<OSteps::FDStep_Type::IMPLICIT>(M,OSteps::make_time(), OSteps::make_context(m)); 
  // cout << M << endl; 

  LinOps::VectorXD vec = LinOps::make_Discretization(m, [](double x, double y){ return x;}); 

  // bcs.VecBeforeStep<OSteps::FDStep_Type::IMPLICIT>(vec.values(), OSteps::make_time(), OSteps::make_context(m)); 
  bcs.VecAfterStep<OSteps::FDStep_Type::EXPLICIT>(vec.values(), OSteps::make_time(), OSteps::make_context(m)); 

  print_mat(m->OneDim_views(vec.values()), "vec"); 
};