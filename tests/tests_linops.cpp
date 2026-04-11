// tests_linops.cpp
//
// linops Test suite with GTest
//
// JAF 12/5/2025

#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include<LinOps/All.hpp>

using namespace linops; 
using namespace linops::internal; 

// Mesh Suite ---------------------------------------- 
TEST(MeshSuite1D, Mesh1DConstructible){
  // simply make a mesh and do nothing with it
  Mesh1D my_mesh; 

  // make a mesh with custom endpoints, # of steps
  Mesh1D my_custom_mesh(-10.0,20.0,31); 
  ASSERT_EQ(my_custom_mesh.size(), 31); 

  // mesh throws errors if # of grid points < 2
  EXPECT_ANY_THROW(Mesh1D(0.0, 1.0, 1));

  // mesh throws errors if x1 > x2 
  EXPECT_ANY_THROW(Mesh1D(1.0, -5.0, 6));
};

TEST(MeshSuite1D, Mesh1DIndexing){
  // make a mesh 
  Mesh1D my_mesh(0.0,10.0,11);  

  ASSERT_EQ(my_mesh[0], 0.0); 
  ASSERT_EQ(my_mesh[10], 10.0); 
  ASSERT_EQ(my_mesh.at(0), 0.0); 
  ASSERT_EQ(my_mesh.at(10), 10.0); 
  
  // operator[] doesn't check size 
  my_mesh[20]; // way out of bounds. but no error
  EXPECT_ANY_THROW(my_mesh.at(30)); 
};

TEST(MeshSuite1D, Mesh1DIterators){
  // simply make a mesh and do nothing with it
  int n_steps = 11;
  Mesh1D my_mesh(0.0,10.0,n_steps);

  // give all iterators as std::vec 
  my_mesh.begin(); 
  my_mesh.end(); 
  my_mesh.cbegin(); 
  my_mesh.cend(); 
  my_mesh.rbegin(); 
  my_mesh.rend(); 
  my_mesh.crbegin(); 
  my_mesh.crend(); 
  
  int count=0; 
  std::vector<double>::iterator it = my_mesh.begin(); 
  while(it!= my_mesh.end()){ count++; it++;}; 
  ASSERT_EQ(count,n_steps); 

  int reverse_count=0; 
  std::vector<double>::reverse_iterator rit = my_mesh.rbegin(); 
  while(rit!= my_mesh.rend()){ reverse_count++; rit++;}; 
  ASSERT_EQ(reverse_count,n_steps); 
};



// Vector Suite ---------------------------------------- 
TEST(VectorSuite1d, Disc1DConstructible)
{
  linops::Vector1D my_vals; 
}

TEST(VectorSuite1d, Disc1DMovable)
{
  int n_steps=11; 
  auto my_mesh = make_Mesh1D(0.0,10.0,n_steps); 
  linops::Vector1D moved_from(my_mesh);
  linops::Vector1D moved_to(std::move(moved_from));  
  // ASSERT_TRUE(moved_from.getMesh1D().expired()); // no longer altering mesh in moved_from 
  ASSERT_EQ(moved_from.values().data(),nullptr); // moved from now has invalid eigen::vectorxd 
  ASSERT_EQ(moved_to.size(), n_steps); 
  ASSERT_EQ(moved_to.getMesh1D(), my_mesh); 
}

TEST(VectorSuite1d, Disc1DSetMesh)
{
  auto my_mesh = make_Mesh1D(); 
  linops::Vector1D my_vals; 
  ASSERT_FALSE(my_vals.getMesh1D()); 
  linops::Vector1D discretization_w_stored_mesh(my_mesh); 
  ASSERT_TRUE(discretization_w_stored_mesh.getMesh1D());
}

TEST(VectorSuite1d, Disc1DSetCosntant)
{
  auto my_mesh = make_Mesh1D(); 
  double val_set = 0.0; 

  linops::Vector1D my_vals = linops::make_Discretization(my_mesh, val_set); 
  
  ASSERT_EQ(my_vals.at(0), val_set); 
  ASSERT_EQ(my_vals.at(my_vals.size()-1), val_set); 
}

TEST(VectorSuite1d, Disc1DSetByCallable)
{
  auto my_lambda = [](const double& x){return x*x;}; 
  struct Callable_t 
  {
    double operator()(const double& x){return x*x*x;}; 
  }; 
  Callable_t my_callable; 

  int n_steps = 101; 
  double left=0, right=100; 
  auto my_mesh = make_Mesh1D(left,right,n_steps); 

  // with lambda 
  linops::Vector1D my_vals = linops::make_Discretization(my_mesh,my_lambda); 

  ASSERT_EQ(my_vals.at(0), my_lambda(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_lambda(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps-1), my_lambda(my_mesh->at(n_steps-1))); 

  my_vals = linops::make_Discretization(my_mesh,my_callable);

  ASSERT_EQ(my_vals.at(0), my_callable(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_callable(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps-1), my_callable(my_mesh->at(n_steps-1))); 

}

TEST(VectorSuite1d, Disc1DIterators)
{
    // simply make a mesh and do nothing with it
  int n_steps = 11;
  auto my_mesh = make_Mesh1D(0.0,10.0,n_steps);

  linops::Vector1D my_vals(my_mesh); 

  // give all iterators as std::vec 
  my_vals.begin(); 
  my_vals.end(); 
  my_vals.cbegin(); 
  my_vals.cend(); 
  // Eigen::VectorXd has no reverse iterators. 
  // my_vals.rbegin(); 
  // my_vals.rend(); 
  // my_vals.crbegin(); 
  // my_vals.crend(); 
  
  int count=0; 
  Eigen::VectorXd::iterator it = my_vals.begin(); 
  while(it!= my_vals.end()){ count++; it++;}; 
  ASSERT_EQ(count,n_steps); 
};



// Linear Operators Suite ---------------------------------------------
// Just the LinOpBase CRTP class. 
TEST(LinearOperatorSuite, IdentityConstructible)
{
  // default construct uses nullptr
  IOp Identity01;
  ASSERT_EQ(Identity01.getMesh1D(),nullptr);

  // construct with ptr arg
  auto my_mesh = make_Mesh1D(); 
  IOp Identity02(my_mesh);
  ASSERT_EQ(Identity02.getMesh1D(), my_mesh);   

  // Check that entries on diag are 1
  int s = my_mesh->size()-1; 
  linops::Matrix M02 = Identity02.asMatrix(); 
  ASSERT_EQ(M02.coeff(0,0),1); 
  ASSERT_EQ(M02.coeff(s,s),1); 
  ASSERT_EQ(M02.coeff(s/2,s/2),1);

  // of diag are zero
  ASSERT_EQ(M02.coeff(0,s),0); 
  ASSERT_EQ(M02.coeff(s,0),0); 
};

TEST(LinearOperatorSuite, RandOpConstructible)
{
  // default construct uses nullptr
  RandOp Rand01;
  ASSERT_EQ(Rand01.getMesh1D(),nullptr);

  // construct with ptr arg
  auto my_mesh = make_Mesh1D(); 
  RandOp Rand02(my_mesh);
  ASSERT_EQ(Rand02.getMesh1D(), my_mesh);   
};

TEST(LinearOperatorSuite, RandOp_asMat)
{
  // construct with ptr arg
  auto my_mesh = make_Mesh1D(); 
  RandOp Rand01(my_mesh);

  ASSERT_EQ(Rand01.getMesh1D(), my_mesh);   

  Eigen::MatrixXd result = RandOp().asMatrix(); 
};

TEST(LinearOperatorSuite, RandOpApply)
{
  // setup mesh + discretization 
  auto my_mesh = make_Mesh1D(); 
  auto func = [](double x){return x*x;}; // x^2 

  linops::Vector1D my_vals = linops::make_Discretization(my_mesh, func);

  // get a random linear operator 
  RandOp Rand01(my_mesh);
  // get its underlying matrix representation 
  linops::Matrix matrix_rep = Rand01.asMatrix();

  // make sure .apply() gives the same as A*v 
  linops::Vector1D apply_method_result = Rand01.apply(my_vals); 
  Eigen::VectorXd manual_linalg_result = matrix_rep * my_vals.values(); 
  for(std::size_t i=0; i<apply_method_result.size(); i++){
    ASSERT_EQ(apply_method_result[i], manual_linalg_result[i]); 
  }
};

// Using LinOpExpr. 
TEST(LinearOperatorSuite, BasicAddition)
{
  auto my_mesh = make_Mesh1D(); 
  // construct without ptr arg
  auto make_I_rval = [my_mesh](){return IOp(my_mesh);}; 
  IOp I_lval(my_mesh);

  auto sum01 = I_lval + I_lval;  
  auto sum02 = I_lval + make_I_rval();
  auto sum03 = make_I_rval() + I_lval;
  auto sum04 = make_I_rval() + make_I_rval(); 

  // lambda to check 4 corners + middle of a matrix expr
  auto check_lambda = [s = my_mesh->size()-1](const auto& expr) -> void
  {
    linops::Matrix Mat = expr.asMatrix(); 
    // Check that entries on diag are 2
    ASSERT_EQ(Mat.coeff(0,0),2); 
    ASSERT_EQ(Mat.coeff(s,s),2); 
    ASSERT_EQ(Mat.coeff(s/2,s/2),2);

    // of diag are zero
    ASSERT_EQ(Mat.coeff(0,s),0); 
    ASSERT_EQ(Mat.coeff(s,0),0); 
  };

  check_lambda(sum01); 
  check_lambda(sum02); 
  check_lambda(sum03); 
  check_lambda(sum04); 

  // eigen sparse matrices have no operator==() .... 
  auto to_dense = [](auto expr) -> Eigen::MatrixXd {
    Eigen::MatrixXd result;
    result = expr.asMatrix();
    return result; 
  };

  // compar expr 1=2, 1=3, 1=4
  ASSERT_EQ(to_dense(sum01),to_dense(sum02));
  ASSERT_EQ(to_dense(sum01),to_dense(sum03));
  ASSERT_EQ(to_dense(sum01),to_dense(sum04));
};

TEST(LinearOperatorSuite, ScalarMultiplication)
{
  auto my_mesh = make_Mesh1D(); 
  RandOp L1(my_mesh);

  // test multiplication with rvals
  2.0* RandOp(my_mesh);

  // multiply by scalar 
  double c = 3.0; 
  auto Expr = c * L1; 

  // evalues the expression object to a dense matrix
  auto to_dense = [](auto& expr) -> Eigen::MatrixXd {
    Eigen::MatrixXd result;
    result = expr.asMatrix().eval();
    return result; 
  };

  // compar expr 1=2, 1=3, 1=4
  ASSERT_EQ(to_dense(Expr), (c*to_dense(L1)).eval());
};

TEST(LinearOperatorSuite, Composition)
{
  // get a mesh
  auto my_mesh = make_Mesh1D(); 

  // vector of [1,1,...,1]
  linops::Vector1D my_vals = linops::make_Discretization(my_mesh, 1.0); 

  // construct without ptr arg
  RandOp L1(my_mesh), L2(my_mesh);

  // composition L1( L2(.) ) 
  auto Expr = L1.compose(L2); 

  // // PRINT ---------------------------
  // using std::cout, std::endl; 
  // cout << "my_vals: " << my_vals.values().transpose() << endl; 
  // cout << "L1 -------" <<endl << L1.asMatrix() << endl; 
  // cout << "L2 -------" <<endl << L2.asMatrix() << endl; 
  // // cout << "expr-----" << endl << Expr.asMatrix() << endl; 
  // // PRINT ---------------------------

  // get underlying Eigen::VectorXd results of .apply() 
  auto expression_result = Expr.apply(my_vals); 
  // composition by hand 
  auto manual_result = L1.apply(L2.apply(my_vals)).values(); 

  // size is ==
  ASSERT_EQ(expression_result.size(), manual_result.size()); 
  // assert each value is == 
  for(int i=0; i< expression_result.size(); i++){
    ASSERT_NEAR(expression_result[i], manual_result[i], 1e-4);
  }
};

TEST(LinearOperatorSuite, ExpressionChaining)
{
  auto my_mesh = make_Mesh1D(); 
  // just a messy expression 
  // auto my_expr = (2.0*(2.0*(2.0*(2.0*IOp())))).compose(50*IOp(my_mesh) + RandOp() + IOp() - RandOp(my_mesh).compose(IOp(my_mesh)));
  // not all lhs/rhs had a mesh in expression construction 
  // my_expr.setMesh1D(my_mesh); 
  // we should be able to make into eigen Matrix no matter what
  // Eigen::MatrixXd resulting_mat = my_expr.asMatrix(); 

  // building the same expression step by step 
  auto tmp1 = 2.0*IOp();
  auto tmp2 = 2.0*tmp1; 
  auto tmp3 = 2.0*tmp2; 
  auto tmp4 = 2.0*tmp3;
  auto my_expr = tmp4; 

  auto rhs1 = 50*IOp(my_mesh); 
  auto rhs2 = rhs1 + RandOp(); 
  auto rhs3 = rhs2 + IOp(); 
  auto rhs4 = RandOp(my_mesh).compose(IOp(my_mesh));
  auto my_expr2 = rhs4; // all temporaries still alive

  // test the expression still has a .apply() method 
  linops::Vector1D disc = linops::make_Discretization(my_mesh, 1.0); 

  linops::Matrix M1 = my_expr.asMatrix(); 
  my_expr2.setMesh1D(my_mesh);
  my_expr2.apply(disc);  

  linops::Matrix M2 = my_expr2.asMatrix(); 
  my_expr2.setMesh1D(my_mesh);
  my_expr2.apply(disc);  
  my_expr2.setMesh1D(my_mesh); 
}

TEST(LinearOperatorSuite, Method_set_mesh_ExprHooking)
{
  // Calling setMesh1D on an expression E = L1 + L2 
  // should pass the mesh to L1.setMesh1D() and L2.setMesh1D() 
   
  // construct without mesh ptrs 
  IOp I_lval;
  auto Expr = I_lval + IOp();

  // make mesh and give it to expression
  auto my_mesh = make_Mesh1D();
  Expr.setMesh1D(my_mesh); 

  // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
  ASSERT_EQ(I_lval.getMesh1D(), my_mesh);
  ASSERT_EQ(Expr.getRhs().getMesh1D(), my_mesh);

  // test again for scalar multiply  // construct without mesh ptrs 
  IOp I2_lval;
  double c=2.0; 
  auto Expr2 = 2.0* I2_lval;
  auto Expr3 = c*IOp();
  Expr2.setMesh1D(my_mesh); 
  Expr3.setMesh1D(my_mesh); 
  // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
  ASSERT_EQ(I2_lval.getMesh1D(), my_mesh);
  ASSERT_EQ(Expr2.getRhs().getMesh1D(), my_mesh);
  ASSERT_EQ(Expr3.getRhs().getMesh1D(), my_mesh);

  // test again for composition 
  RandOp I3_lval; 
  IOp I4_lval;
  auto Expr4 = I3_lval.compose(I3_lval); 
  auto Expr5 = Expr4.compose(IOp());
  // Expr4.setMesh1D(my_mesh); 
  Expr5.setMesh1D(my_mesh); 
  // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
  ASSERT_EQ(I3_lval.getMesh1D(), my_mesh);
  // ASSERT_EQ(Expr4.Rhs().getMesh1D(), my_mesh);
  ASSERT_EQ(Expr5.getRhs().getMesh1D(), my_mesh);
}

/* TEST(LinearOperatorSuite, LinOpTraits)
{
  struct foo {
    Discretization1D apply(const Discretization1D& d) const { return d; };
  };        
  struct bar{};

  // given a potential crtp mixin. see if it has .apply(discretization1d) -> const discretization1d method
  ASSERT_TRUE(traits::has_apply<foo>::value); 
  ASSERT_FALSE(traits::has_apply<bar>::value); 

  // see if a given type is a derived from the LinOpBase<> CRTP class
  ASSERT_TRUE(traits::is_linop_crtp<RandOp>::value);
  ASSERT_FALSE(traits::is_linop_crtp<int>::value); 
}
*/ 

// testing out setTime() hooking
TEST(FdmPluginSuite, Method_SetTime_Hooking)
{
  // make a LinOP
  linops::IOp I; 

  auto expr01 = I+I;
  auto expr02 = I.compose(I); 
  auto expr03 = 3.0*I; 
  auto expr04 = I.compose(linops::IOp()); 

  // set time, make sure it is equal 
  I.setTime(1.0);
  ASSERT_EQ(1.0,I.getTime()); 

  // set time of expression, make sure it propagates to t 
  expr01.setTime(2.0);
  ASSERT_EQ(2.0,I.getTime()); 

  // set time of expression, make sure it propagates to t 
  expr02.setTime(3.0);
  ASSERT_EQ(3.0,I.getTime()); 

  // set time of expression, make sure it propagates to t 
  expr03.setTime(4.0);
  ASSERT_EQ(4.0,I.getTime()); 

  // set time of expression, make sure it propagates to t 
  expr04.setTime(5.0);
  ASSERT_EQ(5.0,I.getTime()); 


  // make sure LHS is given priority 
  linops::IOp I1, I2; 
  I1.setTime(1.0); 
  I2.setTime(2.0); 
  auto expr = I1+I2; 

  ASSERT_EQ(I1.getTime(), expr.getTime()); 
  ASSERT_TRUE(I2.getTime() != expr.getTime()); 
}; 

// testing NthDerivOp is constructible 
TEST(NthDerivOpSuite, NthDerivOpConstructible)
{
  NthDerivOp<1> my_deriv; 
  NthDerivOp<2> order_2; 
  // NthDerivOp<2> from_mesh(nullptr); 
}

// testing setMesh1D() completes with no errors. 
TEST(NthDerivOpSuite, Method_set_mesh_completing)
{
  auto my_mesh_01 = linops::make_Mesh1D(0.0,10.0,11);
  auto my_mesh_02 = linops::make_Mesh1D(0.0,10.0,101);
  auto my_mesh_03 = linops::make_Mesh1D(0.0,10.0,1001);
  auto my_mesh_04 = linops::make_Mesh1D(0.0,10.0,10001);

  auto D1 = NthDerivOp<1>{}; 
  auto D2 = NthDerivOp<2>{}; 
  auto D3 = NthDerivOp<3>{}; 
  auto D4 = NthDerivOp<4>{}; 

  auto test_mesh_lam = [&](auto mesh_ptr){
    D1.setMesh1D(mesh_ptr); 
    D1.asMatrix(); 
    D2.setMesh1D(mesh_ptr); 
    D2.asMatrix(); 
    D3.setMesh1D(mesh_ptr); 
    D3.asMatrix(); 
    D4.setMesh1D(mesh_ptr);
    D1.asMatrix(); 
  };

  test_mesh_lam(my_mesh_01);
  test_mesh_lam(my_mesh_02);
  test_mesh_lam(my_mesh_03);
  test_mesh_lam(my_mesh_04);
}

/* // testing NthDerivOp custom .compose() method  
TEST(NthDerivOpSuite, NthDerivOpCompose)
{
  using D = linops::NthDerivOp; 
  // constructible with a certain order 
  ASSERT_EQ(D(3).Order(), 3); 

  // orders for derivative 
  std::size_t n=4, m=7; 

  // testing composition splits sums ((f+g)' = f' + g')
  auto add_expr = D(n) + D(m); 
  ASSERT_EQ(D(n).compose(add_expr).getLhs().Order(), n+n); 
  ASSERT_EQ(D(n).compose(add_expr).getRhs().Order(), n+m);
  
  // testing composition parses scalar multiply ((c*f)' = c*f')
  auto mult_expr = 4.0 * D(n); 
  // ASSERT_EQ(D(n).compose(mult_expr).getRhs().Order(), n+n);
}

*/ 

