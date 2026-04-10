// tests_texprs.cpp
//
// GTest suite for TExprs classes 
//
// JAF 1/23/2026

#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include<LinOps/All.hpp>
#include<OutsideSteps/All.hpp> 
#include<TExprs/All.hpp>

// DerivExpressionSuite ---------------------------------------------------------
TEST(DerivExpressionSuite, NthTimeDeriv)
{
  TExprs::NthTimeDeriv order_1(1); 
  TExprs::NthTimeDeriv order_2(2); 
  TExprs::NthTimeDeriv order_3(3);
  
  // testing .Order() method 
  ASSERT_EQ(order_1.Order(), 1); 
  ASSERT_EQ(order_2.Order(), 2); 
  ASSERT_EQ(order_3.Order(), 3); 

  // testing .CoeffAt() method
  std::vector<double> weights = {0,1,2,   3,4,5,    6,7,8,   9,10,11}; 
  std::size_t node_per_row = 3; 
  ASSERT_EQ(order_1.CoeffAt(weights,node_per_row ,0), 3.0); 
  ASSERT_EQ(order_1.CoeffAt(weights,node_per_row ,1), 4.0); 
  ASSERT_EQ(order_1.CoeffAt(weights,node_per_row ,2), 5.0); 

  ASSERT_EQ(order_2.CoeffAt(weights,node_per_row ,0), 6.0); 

  ASSERT_EQ(order_3.CoeffAt(weights,node_per_row ,0), 9.0); 

  // testing .toTuple() method 
  auto tup = order_1.toTuple(); 
  ASSERT_EQ(std::get<0>(tup).Order(), 1); 
}

TEST(DerivExpressionSuite, CoeffMultExpr)
{
  TExprs::NthTimeDeriv order_1(1); 
  TExprs::NthTimeDeriv order_2(2); 

  double coeff_01 = 2.0; 
  double coeff_02 = 5.0; 

  TExprs::internal::CoeffMultExpr<decltype(coeff_01), decltype(order_1)> mult_01(coeff_01, order_1);
  TExprs::internal::CoeffMultExpr<decltype(coeff_02), decltype(order_2)> mult_02(coeff_02, order_2); 
  
  // testing .Order() method 
  ASSERT_EQ(mult_01.Order(), 1); 
  ASSERT_EQ(mult_02.Order(), 2); 


  // testing .CoeffAt() method
  std::vector<double> weights = {0,1,2,   3,4,5,    6,7,8}; 
  std::size_t node_per_row = 3; 
  ASSERT_EQ(mult_01.CoeffAt(weights,node_per_row ,0), coeff_01 * 3.0); 
  ASSERT_EQ(mult_01.CoeffAt(weights,node_per_row ,1), coeff_01 * 4.0); 
  ASSERT_EQ(mult_01.CoeffAt(weights,node_per_row ,2), coeff_01 * 5.0); 

  ASSERT_EQ(mult_02.CoeffAt(weights,node_per_row ,0), coeff_02 * 6.0); 

  // testing .CoeffAt() when LHS is a matrix
  // TODO!!!!!!!!!!!

  // testing .toTuple() method 
  auto tup = mult_01.toTuple(); 
  ASSERT_EQ(std::get<0>(tup).Order(), 1); 
}

TEST(DerivExpressionSuite, SumExpr)
{
  TExprs::NthTimeDeriv order_1(1); 
  TExprs::NthTimeDeriv order_2(2); 

  TExprs::internal::SumExpr sum_01(std::tuple_cat(order_1.toTuple(), order_2.toTuple()), 67); 

  // testing .Order() method 
  ASSERT_EQ(sum_01.Order(), 67); 


  // testing .CoeffAt() method 
  // will not compile: deleted member func 
  // std::vector<double> weights = {0,1,2,   3,4,5,    6,7,8}; 
  // sum_01.CoeffAt(weights,3, 0); 

  // testing .toTuple() method 
  auto tup = sum_01.toTuple(); 
  ASSERT_EQ(std::tuple_size<decltype(tup)>::value, 2);  
}

TEST(DerivExpressionSuite, TimeExpressionOperators)
{

}

// ExecutorSuite ---------------------------------------------------------- 
/* TEST(TimeExecutorSuite, ExecutorConstructible)
{

}

TEST(TimeExecutorSuite, StoredTimesVec)
{

}

TEST(TimeExecutorSuite, WeightsCalc)
{

}

TEST(TimeExecutorSuite, StoredSolsVec)
{

}

TEST(TimeExecutorSuite, BuildRHS)
{

}

*/ 

// GenSolverSuite ------------------------------------------------
/* TEST(GenSolverSuite, GenSolverConstructible)
{

}

TEST(GenSolverSuite, UtilityMethods)
{

}

*/ 

/* TEST(GenSolverSuite, Solve1DimPDE)
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // defining Domain Mesh --------------------------------------
  auto r = 3.14159; // pi 
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_Mesh1D(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = LinOps::make_Mesh1D(0.0, 5.0, 8); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::Discretization1D my_vals;
  auto sin_lam = [](double x){return 2.0 * std::sin(x); }; 
  my_vals.set_init(my_mesh, sin_lam); 

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(1); 

  // building RHS expression -----------------------------------------------------
  using D = LinOps::NthDerivOp;
  auto space_expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCPair>(left,right); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(1, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenSolver s(time_expr, space_expr); 
  auto v1 = s.Calculate(args); 
  auto v2 = s.CalculateImp(args); 
}

*/ 

/* TEST(GenSolverSuite, SolveXDimPDE)
{
  // defining Domain Mesh --------------------------------------
  auto r = 3.14159; // pi  
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_MeshXD(0.0,r,n_gridpoints,2); 
  // mesh in time 
  auto time_mesh = LinOps::make_Mesh1D(0.0, 3.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::DiscretizationXD my_vals;

  auto sin_lambda = [](double x, double y){ return 3.0 * std::sin(x) * std::sin(y); }; 
  my_vals.set_init(my_mesh, sin_lambda); 

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(2); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::DirectionalNthDerivOp;
  auto space_expr = 0.2 * D(2,0) + 0.2 * D(2,1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_neumann(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCListXD>();
  bcs->list.emplace_back(left,right); 
  bcs->list.emplace_back(left,right); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(2, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenSolver s(time_expr, space_expr); 
  auto v1 = s.Calculate(args); 
  auto v2 = s.CalculateImp(args); 
}

*/

// GenInterpSuite ----------------------------------------------------------
/* TEST(GenInterpSuite, GenInterpConstructible)
{

}

*/  

/* TEST(GenInterpSuite, XDimInterpComplete)
{
  // defining Domain Mesh --------------------------------------
  auto r = 3.14159; // pi  
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_MeshXD(0.0,r,n_gridpoints,2); 
  // mesh in time 
  auto time_mesh = LinOps::make_Mesh1D(0.0, 3.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::DiscretizationXD my_vals;

  auto sin_lambda = [](double x, double y){ return 3.0 * std::sin(x) * std::sin(y); }; 
  my_vals.set_init(my_mesh, sin_lambda); 

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(2); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::DirectionalNthDerivOp;
  auto space_expr = 0.2 * D(2,0) + 0.2 * D(2,1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_neumann(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCListXD>();
  bcs->list.emplace_back(left,right); 
  bcs->list.emplace_back(left,right); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(2, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenInterp interp(time_expr, space_expr, args); 

  // t=0, x=1, y=1
  double val_01 = interp.SolAt(0.0,1.0,1.0); 

  // t=1, x=1, y=1
  double val_02 = interp.SolAt(1.0,1.0,1.0); 

  // t=2, x=1, y=1
  double val_03 = interp.SolAt(2.0,1.0,1.0); 
}

*/ 

/* TEST(GenInterpSuite, 1DimInterpComplete)
{
  // defining Domain Mesh --------------------------------------
  auto r = 3.14159; // pi  
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_Mesh1D(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = LinOps::make_Mesh1D(0.0, 3.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::Discretization1D my_vals;

  auto sin_lambda = [](double x){ return 3.0 * std::sin(x); }; 
  my_vals.set_init(my_mesh, sin_lambda); 

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(1); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto space_expr = 0.5 * D(2) + 0.2 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_neumann(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCPair>(left,right);

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>{my_vals.values()}, 
    .time_dep_flag = false 
  }; 

  TExprs::GenInterp interp(time_expr, space_expr, args); 

  // t=0, x=1
  double val_01 = interp.SolAt(0.0,1.0); 

  // t=1, x=1
  double val_02 = interp.SolAt(1.0,1.0); 

  // t=2, x=1
  double val_03 = interp.SolAt(2.0,1.0); 
}

*/ 

