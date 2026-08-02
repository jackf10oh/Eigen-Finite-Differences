// tests_linops.cpp 
// 
// test diffops + coeffs subdirectories 
// through their basic interface, operators,
// and traits
// 
// JAF 7/4/2026 

// #define FORNFDM_CUSTOM_SCALAR float
#include<fornfdm/plugin.hpp>
#include<fornfdm/all.hpp>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

using namespace fornfdm; 

// Node Selector Suite ------------------------------------------------- 
template< template< std::size_t> class TagType01, template< std::size_t> class TagType02>
struct TestStruct
{
  template<std::size_t M, std::size_t N>
  struct CheckTags
  {
    using Tag01 = TagType01<M>;
    using Tag02 = TagType02<N>;
    static_assert(linops::internal::promote_node_selector_tags<Tag01, Tag02>::is_match, "node selector tags don't match!");
    using ResultTag = typename linops::internal::promote_node_selector_tags<Tag01, Tag02>::type; 
    using CheckTag = std::conditional_t<(M>N), Tag01, Tag02>;
    static_assert(std::is_same_v<ResultTag, CheckTag>, "node selector doesn't take maximum!");
  };
};

TEST(NodeSelectorSuite, BasicInterface){
  // centered
  TestStruct<linops::Centered, linops::Centered>::CheckTags<0,1>{};
  TestStruct<linops::Centered, linops::Centered>::CheckTags<4,7>{};
  TestStruct<linops::Centered, linops::Centered>::CheckTags<5,2>{};
  // backward
  TestStruct<linops::Backward, linops::Backward>::CheckTags<0,1>{};
  TestStruct<linops::Backward, linops::Backward>::CheckTags<4,7>{};
  TestStruct<linops::Backward, linops::Backward>::CheckTags<5,2>{};
  // forward
  TestStruct<linops::Forward, linops::Forward>::CheckTags<0,1>{};
  TestStruct<linops::Forward, linops::Forward>::CheckTags<4,7>{};
  TestStruct<linops::Forward, linops::Forward>::CheckTags<5,2>{};
  // should break at compile time 
  // TestStruct<linops::Centered, linops::Forward>::CheckTags<0,1>{};
  // TestStruct<linops::Forward, linops::Backward>::CheckTags<4,7>{};
  // TestStruct<linops::Centered, linops::Backward>::CheckTags<5,2>{};
}

// Linops Suite ------------------------------------------------- 
TEST(LinopsSuite, Constructors){
  using namespace linops;

  NthPartialDeriv<0, 0> linop_00;
  NthPartialDeriv<1, 0> linop_01;
  NthPartialDeriv<2, 0> linop_02;
  NthPartialDeriv<3, 0> linop_03;

  NthPartialDeriv<0, 0, Centered<20>> linop_04;
  NthPartialDeriv<1, 0, Centered<20>> linop_05;
  NthPartialDeriv<2, 0, Centered<20>> linop_06;
  NthPartialDeriv<3, 0, Centered<20>> linop_07;

  NthPartialDeriv<0, 1> linop_08;
  NthPartialDeriv<1, 1> linop_09;
  NthPartialDeriv<2, 1> linop_10;
  NthPartialDeriv<3, 1> linop_11;

  NthPartialDeriv<0, 1, Forward<20>> linop_12;
  NthPartialDeriv<1, 1, Forward<20>> linop_13;
  NthPartialDeriv<2, 1, Forward<20>> linop_14;
  NthPartialDeriv<3, 1, Forward<20>> linop_15;

  NthPartialDeriv<0, 2> linop_16;
  NthPartialDeriv<1, 2> linop_17;
  NthPartialDeriv<2, 2> linop_18;
  NthPartialDeriv<3, 2> linop_19;

  NthPartialDeriv<0, 2, Backward<20>> linop_20;
  NthPartialDeriv<1, 2, Backward<20>> linop_21;
  NthPartialDeriv<2, 2, Backward<20>> linop_22;
  NthPartialDeriv<3, 2, Backward<20>> linop_23;
};

TEST(LinopsSuite, BasicInterface){
  linops::NthPartialDeriv<1, 0> Ux;

  // setter + getter to domain
  auto mesh = make_Mesh(fornfdm::linspaced(21,0.0,20.0));
  Ux.setMesh(mesh);
  std::shared_ptr<const fornfdm::Mesh> copy = Ux.getMesh();
  
  // setter + getter to time
  Ux.setTime(0.0);
  fornfdm::Real t_copy = Ux.getTime();
};

TEST(LinopsSuite, BinaryAddition){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;
  auto xpr = Ux + Uxx;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);
  xpr.setMesh(mesh);

  fornfdm::CSRMatrix result = xpr;
  for(auto i=0; i<n; ++i)
  {
    for(auto j=0; j<n; ++j)
    {
      auto val = Ux.getStencil().coeff(i,j) + Uxx.getStencil().coeff(i,j);
      ASSERT_EQ(result.coeff(i,j), val);
    }
  }
};

TEST(LinopsSuite, BinarySubtraction){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;
  auto xpr = Ux - Uxx;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);
  xpr.setMesh(mesh);

  fornfdm::CSRMatrix result = xpr;
  for(auto i=0; i<n; ++i)
  {
    for(auto j=0; j<n; ++j)
    {
      auto val = Ux.getStencil().coeff(i,j) - Uxx.getStencil().coeff(i,j);
      ASSERT_EQ(result.coeff(i,j), val);
    }
  }
};

TEST(LinopsSuite, UnaryNegation){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;
  auto xpr = -Ux;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);
  xpr.setMesh(mesh);

  fornfdm::CSRMatrix result = xpr;
  for(auto i=0; i<n; ++i)
  {
    for(auto j=0; j<n; ++j)
    {
      auto val = -Ux.getStencil().coeff(i,j);
      ASSERT_EQ(result.coeff(i,j), val);
    }
  }
};

TEST(LinopsSuite, ScalarLeftMultiply){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);

  auto test_lam = [&](fornfdm::Scalar c)
  {
    auto xpr01 = c * Ux;
    auto xpr02 = c * Uxx;
    xpr01.setMesh(mesh);
    xpr02.setMesh(mesh);
    fornfdm::CSRMatrix result01 = xpr01;
    fornfdm::CSRMatrix result02 = xpr02;
    for(auto i=0; i<n; ++i)
    {
      for(auto j=0; j<n; ++j)
      {
        auto val01 = c * Ux.getStencil().coeff(i,j);
        auto val02 = c * Uxx.getStencil().coeff(i,j);
        ASSERT_EQ(result01.coeff(i,j), val01);
        ASSERT_EQ(result02.coeff(i,j), val02);
      }
    }
  };
  test_lam(5.0);
  test_lam(8.0);
  test_lam(13.0);
  test_lam(21.0);
};

TEST(LinopsSuite, AutonCoeffProduct1D){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);

  auto test_lam = [&](auto callable)
  {
    linops::AutonomousCoeff c(callable);
    auto xpr01 = c * Ux;
    auto xpr02 = c * Uxx;
    xpr01.setMesh(mesh);
    xpr02.setMesh(mesh);
    fornfdm::CSRMatrix result01 = xpr01;
    fornfdm::CSRMatrix result02 = xpr02;
    for(auto i=0; i<n; ++i)
    {
      fornfdm::Scalar eval = callable(mesh->getAxis(0)[i]);
      for(auto j=0; j<n; ++j)
      {
        auto val01 = eval * Ux.getStencil().coeff(i,j);
        auto val02 = eval * Uxx.getStencil().coeff(i,j);
        ASSERT_EQ(result01.coeff(i,j), val01);
        ASSERT_EQ(result02.coeff(i,j), val02);
      }
    }
  };
  test_lam([](fornfdm::Scalar x){ return x; });
  test_lam([](fornfdm::Scalar x){ return x * x - 9.0; });
  test_lam([](fornfdm::Scalar x){ return x*x*x + 2.0*x*x - 7.5*x + 4.3; });
  test_lam([](fornfdm::Scalar x){ return std::sin(x); });
};

TEST(LinopsSuite, TimeDepProduct1D){
  // using 4 nodes so that expression doesn't change # of nodes
  linops::NthPartialDeriv<1, 0, linops::Forward<4>> Ux;
  linops::NthPartialDeriv<2, 0, linops::Forward<4>> Uxx;

  const std::size_t n = 21;
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(n,0.0,n-1));
  Ux.setMesh(mesh);
  Uxx.setMesh(mesh);

  auto test_lam = [&](auto callable)
  {
    linops::TimeDepCoeff c(callable);
    auto xpr01 = c * Ux;
    auto xpr02 = c * Uxx;
    xpr01.setMesh(mesh);
    xpr02.setMesh(mesh);
    for(const auto& t : {0.0, 3.0, 6.0, 9.0})
    {
      xpr01.setTime(t);
      xpr02.setTime(t);
      fornfdm::CSRMatrix result01 = xpr01;
      fornfdm::CSRMatrix result02 = xpr02;
      for(auto i=0; i<n; ++i)
      {
        fornfdm::Scalar eval = callable(t, mesh->getAxis(0)[i]);
        for(auto j=0; j<n; ++j)
        {
          auto val01 = eval * Ux.getStencil().coeff(i,j);
          auto val02 = eval * Uxx.getStencil().coeff(i,j);
          ASSERT_EQ(result01.coeff(i,j), val01);
          ASSERT_EQ(result02.coeff(i,j), val02);
        }
      }
    }
  };
  test_lam([](fornfdm::Real t, fornfdm::Scalar x){ return t + x; });
  test_lam([](fornfdm::Real t, fornfdm::Scalar x){ return t + t * x + x * x - 9.0; });
  test_lam([](fornfdm::Real t, fornfdm::Scalar x){ return t*x*x + 2.0*t*x - 7.5*x + 4.3; });
  test_lam([](fornfdm::Real t, fornfdm::Scalar x){ return x*std::sin(t); });
};

TEST(LinopsSuite, EvalTimeMatchesSetTime1D){
  linops::TimeDepCoeff c = [](fornfdm::Real t, fornfdm::Scalar x){ return t + x*x; };
  auto xpr = c * linops::NthPartialDeriv<1,0>{} + linops::NthPartialDeriv<2,0>{}; 

  auto test_lam = [&](int nrows, fornfdm::Real t)
  {
    auto mesh = make_Mesh(fornfdm::linspaced(nrows, 0.0, nrows-1),1);
    xpr.setMesh(mesh); 
    fornfdm::CSRMatrix eval_stencil = xpr.evalTime(t); 
    
    xpr.setTime(t); 
    fornfdm::CSRMatrix set_stencil = xpr; 

    for(int i=0; i<nrows; ++i)
    {
      for(int j=0; j<nrows; ++j)
      {
        ASSERT_EQ(eval_stencil.coeff(i,j), set_stencil.coeff(i,j));
      }
    }
  };
  test_lam(20, 0.0); 
  test_lam(20, 1.0); 
  test_lam(20, 4.0);
  test_lam(100, 0.0); 
  test_lam(100, 1.0); 
  test_lam(100, 4.0);
}

TEST(LinopsSuite, EvalTimeMatchesSetTime2D){

  linops::AutonomousCoeff a = [](fornfdm::Scalar x, fornfdm::Scalar y){ return x * y + y;}; 
  auto dir_00 = a * linops::NthPartialDeriv<2,0>{};

  linops::TimeDepCoeff c = [](fornfdm::Real t, fornfdm::Scalar x){ return t + x*x; };
  auto dir_01 = c * linops::NthPartialDeriv<1,1>{} + linops::NthPartialDeriv<2,1>{}; 

  auto xpr = dir_01 - dir_00;

  auto test_lam = [&](int n, int m, fornfdm::Real t)
  {
    auto mesh = make_Mesh(fornfdm::linspaced(n, 0.0, n-1),fornfdm::linspaced(m, 0.0, m-1));
    
    xpr.setMesh(mesh); 
    fornfdm::CSRMatrix eval_stencil = xpr.evalTime(t); 
    
    xpr.setTime(t); 
    fornfdm::CSRMatrix set_stencil = xpr; 

    for(int i=0; i < n*m; ++i)
    {
      for(int j=0; j < n*m; ++j)
      {
        ASSERT_EQ(eval_stencil.coeff(i,j), set_stencil.coeff(i,j));
      }
    }
  };
  test_lam(11, 11, 0.0);
  test_lam(11, 11, 0.0);
  test_lam(11, 11, 0.0); 

  test_lam(11, 21, 0.0);
  test_lam(11, 21, 0.0);
  test_lam(11, 21, 0.0); 

  test_lam(21, 21, 0.0);
  test_lam(21, 21, 0.0);
  test_lam(21, 21, 0.0); 
  
}

// TODO higher dimension tests