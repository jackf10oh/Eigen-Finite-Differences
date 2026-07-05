// tests_fornberg.cpp 
// 
// tests fornberg algorithm against some 
// known finite difference weights
// 
// JAF 5/18/2026 

#include<fornfdm/Mesh.hpp>
#include<fornfdm/utilities/fornberg.hpp>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

// Solver Suite ------------------------------------------------- 
TEST(FornbergSuite, BasicNodes){
  fornfdm::Vector nodes = fornfdm::linspaced(11,0.0,10.0); 
  std::array<fornfdm::Scalar, 50> result; 

  auto validate_weights = [&](const fornfdm::Scalar* v, std::size_t n, bool debug=false)
  {
    for(auto i=0; i<n; ++i)
    {
      if(debug) std::cout << i << ", "; 
      ASSERT_NEAR(result[i],v[i],1e-9);
    }
    if(debug) std::cout << "\n";
  };

  // forward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[0], 1, result.begin()); 
  fornfdm::Scalar forward_first[4] = {1.0, 0.0, -1.0, 1.0}; 
  validate_weights(forward_first, 4); 

  // backward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[1], 1, result.begin()); 
  fornfdm::Scalar backward_first[4] = {0.0, 1.0, -1.0, 1.0}; 
  validate_weights(backward_first, 4); 

  // centered. 3 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+3, nodes[1], 2, result.begin()); 
  fornfdm::Scalar centered_second[9] = {0.0, 1.0, 0.0, -0.5, 0.0, 0.5, 1, -2, 1}; 
  validate_weights(centered_second, 9); 

  // centered. 5 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+5, nodes[2], 2, result.begin()); 
  fornfdm::Scalar centered_second_05[15] = {
    0, 0, 1, 0, 0,
    1.0/12, -2.0/3, 0, 2.0/3, -1.0/12, 
    -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12, 
  }; 
  validate_weights(centered_second_05, 15); 

  // forward. 4 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+4, nodes[0], 2, result.begin()); 
  fornfdm::Scalar forward_second[12] = {
    1, 0, 0, 0,
    -11.0/6, 3, -3.0/2, 1.0/3, 
    2, -5, 4, -1
  }; 
  validate_weights(forward_second, 12); 

  // forward. 9 nodes. 3rd order -> 45 entries
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+9, nodes[0], 3, result.begin()); 
  fornfdm::Scalar forward_third_09[45] = {
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    -761.0/280, 8, -14, 56.0/3, -35.0/2, 56.0/5, -14.0/3, 8.0/7, -1.0/8, 
    29531.0/5040, -962.0/35, 621.0/10, -4006.0/45, 691.0/8, -282.0/5, 2143.0/90, -206.0/35, 363.0/560, 
    -801.0/80, 349.0/6, -18353.0/120, 2391.0/10, -1457.0/6, 4891.0/30, -561.0/8, 527.0/30, -469.0/240
  }; 
  validate_weights(forward_third_09, 45); 
};

TEST(FornbergSuite, ShiftedNodes01){
  fornfdm::Vector nodes = fornfdm::linspaced(11,1.0,11.0); 
  std::array<fornfdm::Scalar, 50> result; 

  auto validate_weights = [&](const fornfdm::Scalar* v, std::size_t n, bool debug=false)
  {
    for(auto i=0; i<n; ++i)
    {
      if(debug) std::cout << i << ", "; 
      ASSERT_NEAR(result[i],v[i],1e-9);
    }
    if(debug) std::cout << "\n";
  };

  // forward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[0], 1, result.begin()); 
  fornfdm::Scalar forward_first[4] = {1.0, 0.0, -1.0, 1.0}; 
  validate_weights(forward_first, 4); 

  // backward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[1], 1, result.begin()); 
  fornfdm::Scalar backward_first[4] = {0.0, 1.0, -1.0, 1.0}; 
  validate_weights(backward_first, 4); 

  // centered. 3 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+3, nodes[1], 2, result.begin()); 
  fornfdm::Scalar centered_second[9] = {0.0, 1.0, 0.0, -0.5, 0.0, 0.5, 1, -2, 1}; 
  validate_weights(centered_second, 9); 

  // centered. 5 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+5, nodes[2], 2, result.begin()); 
  fornfdm::Scalar centered_second_05[15] = {
    0, 0, 1, 0, 0,
    1.0/12, -2.0/3, 0, 2.0/3, -1.0/12, 
    -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12, 
  }; 
  validate_weights(centered_second_05, 15); 

  // forward. 4 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+4, nodes[0], 2, result.begin()); 
  fornfdm::Scalar forward_second[12] = {
    1, 0, 0, 0,
    -11.0/6, 3, -3.0/2, 1.0/3, 
    2, -5, 4, -1
  }; 
  validate_weights(forward_second, 12); 

  // forward. 9 nodes. 3rd order -> 45 entries
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+9, nodes[0], 3, result.begin()); 
  fornfdm::Scalar forward_third_09[45] = {
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    -761.0/280, 8, -14, 56.0/3, -35.0/2, 56.0/5, -14.0/3, 8.0/7, -1.0/8, 
    29531.0/5040, -962.0/35, 621.0/10, -4006.0/45, 691.0/8, -282.0/5, 2143.0/90, -206.0/35, 363.0/560, 
    -801.0/80, 349.0/6, -18353.0/120, 2391.0/10, -1457.0/6, 4891.0/30, -561.0/8, 527.0/30, -469.0/240
  }; 
  validate_weights(forward_third_09, 45); 
};

TEST(FornbergSuite, ShiftedNodes02){
  fornfdm::Vector nodes = fornfdm::linspaced(11,-4.0,6.0); 
  std::array<fornfdm::Scalar, 50> result; 

  auto validate_weights = [&](const fornfdm::Scalar* v, std::size_t n, bool debug=false)
  {
    for(auto i=0; i<n; ++i)
    {
      if(debug) std::cout << i << ", "; 
      ASSERT_NEAR(result[i],v[i],1e-9);
    }
    if(debug) std::cout << "\n";
  };

  // forward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[0], 1, result.begin()); 
  fornfdm::Scalar forward_first[4] = {1.0, 0.0, -1.0, 1.0}; 
  validate_weights(forward_first, 4); 

  // backward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[1], 1, result.begin()); 
  fornfdm::Scalar backward_first[4] = {0.0, 1.0, -1.0, 1.0}; 
  validate_weights(backward_first, 4); 

  // centered. 3 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+3, nodes[1], 2, result.begin()); 
  fornfdm::Scalar centered_second[9] = {0.0, 1.0, 0.0, -0.5, 0.0, 0.5, 1, -2, 1}; 
  validate_weights(centered_second, 9); 

  // centered. 5 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+5, nodes[2], 2, result.begin()); 
  fornfdm::Scalar centered_second_05[15] = {
    0, 0, 1, 0, 0,
    1.0/12, -2.0/3, 0, 2.0/3, -1.0/12, 
    -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12, 
  }; 
  validate_weights(centered_second_05, 15); 

  // forward. 4 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+4, nodes[0], 2, result.begin()); 
  fornfdm::Scalar forward_second[12] = {
    1, 0, 0, 0,
    -11.0/6, 3, -3.0/2, 1.0/3, 
    2, -5, 4, -1
  }; 
  validate_weights(forward_second, 12); 

  // forward. 9 nodes. 3rd order -> 45 entries
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+9, nodes[0], 3, result.begin()); 
  fornfdm::Scalar forward_third_09[45] = {
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    -761.0/280, 8, -14, 56.0/3, -35.0/2, 56.0/5, -14.0/3, 8.0/7, -1.0/8, 
    29531.0/5040, -962.0/35, 621.0/10, -4006.0/45, 691.0/8, -282.0/5, 2143.0/90, -206.0/35, 363.0/560, 
    -801.0/80, 349.0/6, -18353.0/120, 2391.0/10, -1457.0/6, 4891.0/30, -561.0/8, 527.0/30, -469.0/240
  }; 
  validate_weights(forward_third_09, 45); 
};

TEST(FornbergSuite, ScaledNodes01){
  fornfdm::Vector nodes = fornfdm::linspaced(11,0.0,100.0); 
  std::array<fornfdm::Scalar, 50> result; 

  auto validate_weights = [&](const fornfdm::Scalar* v, std::size_t n_nodes, std::size_t order, fornfdm::Scalar scale, bool debug=false)
  {
    fornfdm::Scalar rolling = 1.0; 
    for(auto o = 0; o<order+1; ++o)
    {
      for(auto i=0; i<n_nodes; ++i)
      {
        if(debug) std::cout << i << ", "; 
        ASSERT_NEAR(result[o * n_nodes + i],v[o * n_nodes + i] / rolling,1e-9);
      }
      if(debug) std::cout << "\n";
      rolling *= scale;
    }
  };

  // forward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[0], 1, result.begin()); 
  fornfdm::Scalar forward_first[4] = {1.0, 0.0, -1.0, 1.0}; 
  validate_weights(forward_first, 2, 1, 10.0); 

  // backward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[1], 1, result.begin()); 
  fornfdm::Scalar backward_first[4] = {0.0, 1.0, -1.0, 1.0}; 
  validate_weights(backward_first, 2, 1, 10.0); 

  // centered. 3 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+3, nodes[1], 2, result.begin()); 
  fornfdm::Scalar centered_second[9] = {0.0, 1.0, 0.0, -0.5, 0.0, 0.5, 1, -2, 1}; 
  validate_weights(centered_second, 3, 2, 10.0); 

  // centered. 5 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+5, nodes[2], 2, result.begin()); 
  fornfdm::Scalar centered_second_05[15] = {
    0, 0, 1, 0, 0,
    1.0/12, -2.0/3, 0, 2.0/3, -1.0/12, 
    -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12, 
  }; 
  validate_weights(centered_second_05, 5, 2, 10.0); 

  // forward. 4 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+4, nodes[0], 2, result.begin()); 
  fornfdm::Scalar forward_second[12] = {
    1, 0, 0, 0,
    -11.0/6, 3, -3.0/2, 1.0/3, 
    2, -5, 4, -1
  }; 
  validate_weights(forward_second, 4, 2, 10.0); 

  // forward. 9 nodes. 3rd order -> 36 entries
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+9, nodes[0], 3, result.begin()); 
  fornfdm::Scalar forward_third_09[36] = {
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    -761.0/280, 8, -14, 56.0/3, -35.0/2, 56.0/5, -14.0/3, 8.0/7, -1.0/8, 
    29531.0/5040, -962.0/35, 621.0/10, -4006.0/45, 691.0/8, -282.0/5, 2143.0/90, -206.0/35, 363.0/560, 
    -801.0/80, 349.0/6, -18353.0/120, 2391.0/10, -1457.0/6, 4891.0/30, -561.0/8, 527.0/30, -469.0/240
  }; 
  validate_weights(forward_third_09, 9, 3, 10.0); 
};

TEST(FornbergSuite, ScaledNodes02){
  fornfdm::Vector nodes = fornfdm::linspaced(11,0.0,1.0); 
  std::array<fornfdm::Scalar, 50> result; 

  auto validate_weights = [&](const fornfdm::Scalar* v, std::size_t n_nodes, std::size_t order, fornfdm::Scalar scale, bool debug=false)
  {
    fornfdm::Scalar rolling = 1.0; 
    for(auto o = 0; o<order+1; ++o)
    {
      for(auto i=0; i<n_nodes; ++i)
      {
        if(debug) std::cout << i << ", "; 
        ASSERT_NEAR(result[o * n_nodes + i],v[o * n_nodes + i] / rolling,1e-9);
      }
      if(debug) std::cout << "\n";
      rolling *= scale;
    }
  };

  // forward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[0], 1, result.begin()); 
  fornfdm::Scalar forward_first[4] = {1.0, 0.0, -1.0, 1.0}; 
  validate_weights(forward_first, 2, 1, 0.1); 

  // backward stencil. 2 nodes 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+2, nodes[1], 1, result.begin()); 
  fornfdm::Scalar backward_first[4] = {0.0, 1.0, -1.0, 1.0}; 
  validate_weights(backward_first, 2, 1, 0.1); 

  // centered. 3 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+3, nodes[1], 2, result.begin()); 
  fornfdm::Scalar centered_second[9] = {0.0, 1.0, 0.0, -0.5, 0.0, 0.5, 1, -2, 1}; 
  validate_weights(centered_second, 3, 2, 0.1); 

  // centered. 5 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+5, nodes[2], 2, result.begin()); 
  fornfdm::Scalar centered_second_05[15] = {
    0, 0, 1, 0, 0,
    1.0/12, -2.0/3, 0, 2.0/3, -1.0/12, 
    -1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12, 
  }; 
  validate_weights(centered_second_05, 5, 2, 0.1); 

  // forward. 4 nodes. 2nd order 
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+4, nodes[0], 2, result.begin()); 
  fornfdm::Scalar forward_second[12] = {
    1, 0, 0, 0,
    -11.0/6, 3, -3.0/2, 1.0/3, 
    2, -5, 4, -1
  }; 
  validate_weights(forward_second, 4, 2, 0.1); 

  // forward. 9 nodes. 3rd order -> 36 entries
  fornfdm::utils::fornberg(nodes.cbegin(), nodes.cbegin()+9, nodes[0], 3, result.begin()); 
  fornfdm::Scalar forward_third_09[36] = {
    1, 0, 0, 0, 0, 0, 0, 0, 0,
    -761.0/280, 8, -14, 56.0/3, -35.0/2, 56.0/5, -14.0/3, 8.0/7, -1.0/8, 
    29531.0/5040, -962.0/35, 621.0/10, -4006.0/45, 691.0/8, -282.0/5, 2143.0/90, -206.0/35, 363.0/560, 
    -801.0/80, 349.0/6, -18353.0/120, 2391.0/10, -1457.0/6, 4891.0/30, -561.0/8, 527.0/30, -469.0/240
  }; 
  validate_weights(forward_third_09, 9, 3, 0.1); 
};