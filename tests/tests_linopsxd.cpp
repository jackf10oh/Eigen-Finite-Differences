// tests_linopsxd.cpp (LinOpsXD)
//
// Gtest suites for LinOpsXD
//
// JAF 1/5/2026 

#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include<FiniteDifference/LinOps/All.hpp>

using namespace fdm; 
using namespace linops; 
using namespace linops::internal; 

// Mesh Suite ==================================================== 
// Test all constructors 
TEST(MeshXDSuite, MeshXDConstructible){
  // simply make a mesh and do nothing with it
  MeshXD my_mesh; 

  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 
  
  // make a mesh with custom endpoints, # of steps, # of dimension 
  MeshXD from_ends_steps_dims(-10.0,20.0,31, 3); 

  // from std::vec of [left,right] and n_steps 
  MeshXD from_ends_list_steps_list01({{0.0,4.0}}, {8}); 
  MeshXD from_ends_list_steps_list02({{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}, {8,16,32}); 

  // from std::vector<> of std::shared_ptr<Mesh1D> 
  // const!  
  auto mesh_1d_01 = fdm::make_Mesh1D(); 
  auto mesh_1d_02 = fdm::make_Mesh1D(-10.0,0.0,21); 
  auto mesh_1d_03 = fdm::make_Mesh1D(0.0,10.0,101);
  auto v = std::vector<fdm::SharedConstMesh1D>({mesh_1d_01,mesh_1d_02,mesh_1d_03}); 
  MeshXD from_vec(v); 

  // copy 
  MeshXD from_copy(from_dims_only); 

  // // mesh throws errors if # of grid points < 2
  EXPECT_ANY_THROW(MeshXD(0.0, 1.0, 1, 4));

  // // mesh throws errors if x1 > x2 
  EXPECT_ANY_THROW(MeshXD(1.0, -5.0, 6, 1));
}; 

// Testing all method that return std::shared_ptr<Mesh1D> 
TEST(MeshXDSuite, MeshXDMesh1DGetters){
  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 
  const MeshXD from_dims_only_const(X); 

  auto resulting_mesh1d = from_dims_only.getMesh1D(0);
  // !!! not safe to index 
  from_dims_only.getMesh1D(X+5); 
  EXPECT_ANY_THROW(resulting_mesh1d = from_dims_only.getMesh1DSafe(X+5));  

 
  resulting_mesh1d = from_dims_only.getMesh1DSafe(0); 
  // safe to index 
  EXPECT_ANY_THROW(from_dims_only.getMesh1DSafe(X+5));  
}; 

// Testing all method that return sizes 
TEST(MeshXDSuite, MeshXDSizesGetters){
  // simply make a mesh and do nothing with it
  MeshXD my_mesh; 
  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 
  // make a mesh with custom endpoints, # of steps, # of dimension 
  std::size_t n_steps=31; 
  MeshXD from_ends_steps_dims(-10.0,20.0,n_steps, 3); 
  // from std::vec of [left,right] and n_steps 
  std::size_t X1{8}, X2{16}, X3{32}; 
  MeshXD from_ends_list_steps_list01({{0.0,4.0}}, {X1}); 
  MeshXD from_ends_list_steps_list02({{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}, {X1,X2,X3});

  // returns # of dimensions 
  ASSERT_EQ(from_dims_only.numDims(), X); 
  
  // returns # of steps in a specific dimension 
  ASSERT_EQ(from_ends_steps_dims.sizeOfDim(0), n_steps);
  ASSERT_EQ(from_ends_steps_dims.sizeOfDim(1), n_steps);
  EXPECT_ANY_THROW(from_ends_steps_dims.sizeOfDim(10)); // 10 > # of numDims 

  // returns product of (# of steps per dim)
  ASSERT_EQ(from_ends_list_steps_list01.sizesProduct(),X1); 
  ASSERT_EQ(from_ends_list_steps_list02.sizesProduct(),X1*X2*X3);
  
  // returns product of (# of steps per dim) for dim [start,end)
  ASSERT_EQ(from_ends_list_steps_list02.sizesMiddleProduct(0,3), X1*X2*X3); 
  ASSERT_EQ(from_ends_list_steps_list02.sizesMiddleProduct(0,2), X1*X2); 
  ASSERT_EQ(from_ends_list_steps_list02.sizesMiddleProduct(1,3), X2*X3); 
  EXPECT_ANY_THROW(from_ends_list_steps_list02.sizesMiddleProduct(0,10)); // 10 > # of numDims 
  EXPECT_ANY_THROW(from_ends_list_steps_list02.sizesMiddleProduct(4,0)); // start > end  
}; 

// VectorXD Suite ============================================= 
// Test all constructors 
TEST(VectorXDSuite, VectorXDConstructible){
  // default constructor leaves m_vals, m_dims, m_mesh_ptr all empty 
  fdm::VectorXD empty_disc; 

  // 1 dim from single size
  std::size_t N=67; 
  fdm::VectorXD from_size_1d(N); 

  // from X dimensional MeshXDPtr_t
  auto meshes_1d = std::make_shared<MeshXD>(std::size_t{1}); 
  auto meshes_3d = std::make_shared<MeshXD>(std::size_t{3}); 
  fdm::VectorXD from_meshes_1d(meshes_1d); 
  fdm::VectorXD from_meshes_3d(meshes_3d); 

  // copy 
  fdm::VectorXD from_other(from_meshes_1d); 

}

// testing move constructors ???????????????

// Testing all method that return sizes. should match original MeshXDPtr_t 
TEST(VectorXDSuite, VectorXDSizesGetters){
  // simply make a mesh and do nothing with it
  auto my_mesh = make_MeshXD(); 
  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  auto from_dims_only = make_MeshXD(X); 
  // make a mesh with custom endpoints, # of steps, # of dimension 
  std::size_t n_steps=31; 
  auto from_ends_steps_dims = make_MeshXD(-10.0,20.0,n_steps, 3); 

  // from std::vec of [left,right] and n_steps 
  std::size_t X1{8}, X2{16}, X3{32}; 
  using bounds_vec_t = std::vector<std::pair<double,double>>; 
  bounds_vec_t bound_vals01{{0.0,4.0}}, bound_vals02{{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}; 
  using nsteps_vec_t = std::vector<std::size_t>; 
  nsteps_vec_t nsteps01{X1}, nsteps02{X1, X2, X3}; 
  auto from_ends_list_steps_list01 = make_MeshXD(bound_vals01, nsteps01); 
  auto from_ends_list_steps_list02 = make_MeshXD(bound_vals02, nsteps02); 

  // size getters should give same result across VectorXD / MeshXDPtr_t 
  // returns # of dimensions 
  ASSERT_EQ(from_dims_only->numDims(), \
            VectorXD(from_dims_only).numDims()); 
  
  // returns # of steps in a specific dimension 
  ASSERT_EQ(from_ends_steps_dims->sizeOfDim(0), \
            fdm::VectorXD(from_ends_steps_dims).sizeOfDim(0)); 
            
  ASSERT_EQ(from_ends_steps_dims->sizeOfDim(1), \
            fdm::VectorXD(from_ends_steps_dims).sizeOfDim(1));  

  EXPECT_ANY_THROW(fdm::VectorXD(from_ends_steps_dims).sizeOfDim(10)); // 10 > # of numDims 

  // // returns product of (# of steps per dim)
  ASSERT_EQ(from_ends_list_steps_list01->sizesProduct(), \
            fdm::VectorXD(from_ends_list_steps_list01).sizesProduct());  

  ASSERT_EQ(from_ends_list_steps_list02->sizesProduct(), \
            fdm::VectorXD(from_ends_list_steps_list02).sizesProduct());  
  
  // returns product of (# of steps per dim) for dim [start,end)
  ASSERT_EQ(from_ends_list_steps_list02->sizesMiddleProduct(0,3), \
            fdm::VectorXD(from_ends_list_steps_list02).sizesMiddleProduct(0,3));  

  ASSERT_EQ(from_ends_list_steps_list02->sizesMiddleProduct(0,2), \
            fdm::VectorXD(from_ends_list_steps_list02).sizesMiddleProduct(0,2));

  ASSERT_EQ(from_ends_list_steps_list02->sizesMiddleProduct(1,3), \
            fdm::VectorXD(from_ends_list_steps_list02).sizesMiddleProduct(1,3));
 
  EXPECT_ANY_THROW(fdm::VectorXD(from_ends_list_steps_list02).sizesMiddleProduct(0,10)); // 10 > # of numDims 
  EXPECT_ANY_THROW(fdm::VectorXD(from_ends_list_steps_list02).sizesMiddleProduct(4,0)); // start > end 
}

// Testing set_init() for constant
TEST(VectorXDSuite, VectorXDSetByConstant){

  auto my_meshes = make_MeshXD(); 
  
  double val_init = 6.7; 
  auto my_disc = fdm::make_Discretization(my_meshes, val_init); 

  for(auto& val : my_disc.values()){
    ASSERT_EQ(val, val_init); 
  }
}; 

// Testing set_init() for callable  
TEST(VectorXDSuite, VectorXDSetByCallable){
  auto lam00 = [](){return 67.0;}; 
  // Euclidean norms in 1D, 2D, 3D 
  auto lam01 = [](double x){return std::sqrt(x*x);}; 
  auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
  auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

  auto my_mesh_1d = make_MeshXD(0.0,10.0,21, 1); 
  auto my_mesh_2d = make_MeshXD(0.0,10.0,21, 2); 
  auto my_mesh_3d = make_MeshXD(0.0,10.0,21, 3); 

  fdm::VectorXD my_disc; 

  // set init 1D case 
  my_disc = fdm::make_Discretization(my_mesh_1d, lam01); 
  for(std::size_t i=0; i<my_disc.sizeOfDim(0); i++){
    double x = my_mesh_1d->getMesh1DSafe(0)->at(i); 
    double lam_val = lam01(x); 
    double disc_val = my_disc.values()[i];  

    ASSERT_NEAR(disc_val, lam_val, 1e-4); 
  }

  // set init 2D case 
  my_disc = fdm::make_Discretization(my_mesh_2d, lam02); 
  for(std::size_t i=0; i<my_disc.sizeOfDim(0); i++){

    double x = my_mesh_2d->getMesh1DSafe(0)->at(i); 

    for(std::size_t j=0; j<my_disc.sizeOfDim(1); j++){

      double y = my_mesh_2d->getMesh1DSafe(1)->at(j); 

      double lam_val = lam02(x,y); 

      std::size_t flat_idx = i + j*my_disc.sizesMiddleProduct(0,1); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 

    }
  }

  // set init 3D case 
  my_disc = fdm::make_Discretization(my_mesh_3d, lam03); 
  for(std::size_t i=0; i<my_disc.sizeOfDim(0); i++){

    double x = my_mesh_3d->getMesh1DSafe(0)->at(i); 

    for(std::size_t j=0; j<my_disc.sizeOfDim(1); j++){

      double y = my_mesh_3d->getMesh1DSafe(1)->at(j); 

      for(std::size_t k=0; k<my_disc.sizeOfDim(1); k++){

      double z = my_mesh_3d->getMesh1DSafe(2)->at(k); 

      double lam_val = lam03(x,y,z); 

      std::size_t flat_idx = i + j*my_disc.sizesMiddleProduct(0,1) + k*my_disc.sizesMiddleProduct(0,2); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 
      }
    }
  }

  // last case: we can use lower dimensional function on higher dimension mesh
  // i.e. lam02 takes 2 numDims but my_mesh3d has 3 numDims
  my_disc = fdm::make_Discretization(my_mesh_3d, lam02); 
  for(std::size_t i=0; i<my_disc.sizeOfDim(0); i++){

    double x = my_mesh_3d->getMesh1DSafe(0)->at(i); 

    for(std::size_t j=0; j<my_disc.sizeOfDim(1); j++){

      double y = my_mesh_3d->getMesh1DSafe(1)->at(j); 

      for(std::size_t k=0; k<my_disc.sizeOfDim(1); k++){

      double z = my_mesh_3d->getMesh1DSafe(2)->at(k); 

      double lam_val = lam02(x,y); 

      std::size_t flat_idx = i + j*my_disc.sizesMiddleProduct(0,1) + k*my_disc.sizesMiddleProduct(0,2); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 
      }
    }
  }
}; 

// testing move assignment from Eigen::VectorXd ?????????????????

// Linear Operators XD Suite ====================================================

// testing AutonomousCoeff class
TEST(XDOperatorsTestSuite, AutonomousCoeffTest)
{
  // make a mesh 
  auto my_mesh = fdm::make_Mesh1D(0.0,4.0,5); 

  // some lambdas to test out 
  auto lam01 = [](double x){return x*x;}; // x^2
  auto lam02 = [](double x){return std::sin(x);}; // sin(x)
  auto lam03 = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // some polynomial in x 
  
  
  // take 1 coeff_op and 1 lambda. 
  auto check_lam = [&](const auto& coeff_op, auto func){
    Matrix A(coeff_op.asMatrix()); 
    // check each value of diag = func(t,x) or func(t) 
    for(std::size_t i=0; i<A.rows(); i++){
      ASSERT_EQ(    
        A.coeff(i,i), 
        func(my_mesh->at(i))    
      ); 
    }
    // A is NxN where n = mesh size
    ASSERT_EQ(A.rows(),my_mesh->size()); 
    ASSERT_EQ(A.cols(),my_mesh->size());
    
    // of diag is zero 
    ASSERT_EQ(A.coeff(0, A.cols()-1), 0.0); 
    ASSERT_EQ(A.coeff(A.rows()-1, 0), 0.0); 
  }; 

  // create AutonomoousCoeff
  AutonomousCoeff coeff01(lam01,my_mesh); 
  // // check they have the same values
  check_lam(coeff01, lam01);
  

  // create AutonomoousCoeff
  AutonomousCoeff coeff02(lam02,my_mesh); 
  // // check they have the same values
  check_lam(coeff02, lam02);
  
  // create AutonomoousCoeff
  AutonomousCoeff coeff03(lam03,my_mesh); 
  // // check they have the same values
  check_lam(coeff03, lam03);
  
}

// testing TimeDepCoeff class
TEST(CoeffOpTestSuite, TimeDepCoeffTest)
{
  // make a mesh 
  auto my_mesh = fdm::make_Mesh1D(0.0,4.0,5); 

  // some lambdas to test out 
  auto lam01 = [](double t){return t*t;}; // t^2
  auto lam02 = [](double t,double x){return t * std::sin(x);}; // t*sin(x)
  auto lam03 = [](double t,double x){return t*t+t*x+3.0*t*t*t+2*x*x*x-5*x*x+3*x-1;}; // some polynomial in t,x
  
  // take 1 coeff_op and 1 lambda. check each value of diag = func(t,x) or func(t) 
  auto check_lam = [&](const auto& coeff_op, auto func){
    Matrix A(coeff_op.asMatrix()); 
    for(std::size_t i=0; i<A.rows(); i++){
      if constexpr(linops::traits::callable_traits<decltype(func)>::num_args==2){
        ASSERT_EQ(     
          A.coeff(i,i), 
          func(coeff_op.getTime(), my_mesh->cbegin()[i])         
        ); 
      } 
      if constexpr(linops::traits::callable_traits<decltype(func)>::num_args==1){
        ASSERT_EQ(    
          A.coeff(i,i), 
          func(coeff_op.getTime())   
        ); 
      } 
    }
  }; 

  // make coeff
  double t = 3.0; 
  linops::TimeDepCoeff coeff01(lam01,my_mesh); 
  coeff01.setTime(t); 

  // set to new functions / discretization 
  linops::TimeDepCoeff coeff02(lam02,my_mesh); 
  coeff02.setTime(t); 

  // check they have the same values 
  check_lam(coeff02, lam02); 
  
  // set to new functions / discretization 
  linops::TimeDepCoeff coeff03(lam03,my_mesh); 
  coeff03.setTime(t); 

  // check they have the same values 
  check_lam(coeff03, lam03); 
  
} 


