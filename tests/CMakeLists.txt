add_executable(tests_solvers
  tests_solvers.cpp
)

target_link_libraries(tests_solvers
  PRIVATE fornfdm
  PRIVATE GTest::gtest 
  PRIVATE GTest::gtest_main 
  PRIVATE pthread # required by GTest
)

add_executable(tests_fornberg
  tests_fornberg.cpp
)

target_link_libraries(tests_fornberg
  PRIVATE fornfdm
  PRIVATE GTest::gtest 
  PRIVATE GTest::gtest_main 
  PRIVATE pthread # required by GTest
)

include(GoogleTest)
gtest_discover_tests(tests_solvers)
gtest_discover_tests(tests_fornberg)
