add_executable(tests_linops
  tests_linops.cpp
)

target_link_libraries(tests_linops
  PRIVATE FiniteDifference
  PRIVATE GTest::gtest 
  PRIVATE GTest::gtest_main 
  PRIVATE pthread # required by GTest
)

add_executable(tests_linopsxd 
  tests_linopsxd.cpp
)

target_link_libraries(tests_linopsxd 
  PRIVATE FiniteDifference
  PRIVATE GTest::gtest 
  PRIVATE GTest::gtest_main 
  PRIVATE pthread # required by GTest
)

# add_executable(tests_texprs 
#   tests_texprs.cpp
# )

# target_link_libraries(tests_texprs 
#   PRIVATE TExprs
#   PRIVATE GTest::gtest 
#   PRIVATE GTest::gtest_main 
#   PRIVATE pthread # required by GTest
# )

include(GoogleTest)
gtest_discover_tests(tests_linops)
gtest_discover_tests(tests_linopsxd)
# gtest_discover_tests(tests_texprs)
