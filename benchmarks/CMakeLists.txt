add_executable(benchmarks
  benchmarks.cpp
) 

target_link_libraries(benchmarks
  PRIVATE fornfdm 
  PRIVATE benchmark::benchmark
  PRIVATE GTest::gtest 
  PRIVATE GTest::gtest_main 
  PRIVATE pthread # required by GTest
)

if(USE_OPENMP)
  target_link_libraries(benchmarks
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()