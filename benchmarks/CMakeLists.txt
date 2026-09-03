add_executable(benchmarks
  benchmarks.cpp
) 

target_compile_definitions(benchmarks PRIVATE NDEBUG)

target_link_libraries(benchmarks
  PRIVATE fornfdm 
  PRIVATE benchmark::benchmark
  PRIVATE GTest::gmock_main
  # PRIVATE pthread # required by GTest # unnecessary.
)

if(USE_OPENMP)
  target_link_libraries(benchmarks
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()

add_executable(bar_chart
  bar_chart.cpp
) 

target_compile_definitions(bar_chart PRIVATE NDEBUG)

target_link_libraries(bar_chart
  PRIVATE fornfdm 
  PRIVATE benchmark::benchmark
  PRIVATE GTest::gmock_main
  # PRIVATE pthread # required by GTest # unnecessary.
)

if(USE_OPENMP)
  target_link_libraries(bar_chart
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()