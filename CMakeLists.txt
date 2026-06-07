cmake_minimum_required(VERSION 3.16)
project(fornfdm LANGUAGES CXX)
add_library(fornfdm INTERFACE)

# set cxx standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# options 
option(USE_OPENMP "Enable OpenMP Parallelization" OFF) 

# adding OpenMP is enable 
if(USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP found successfully!")
    target_link_libraries(fornfdm INTERFACE OpenMP::OpenMP_CXX)
  endif()
else()
  message(STATUS "Cmake args disabling OpenMP...")
endif()

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# header only library ===============================
target_include_directories(fornfdm 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(fornfdm
  INTERFACE Eigen3::Eigen
)

# adding google test
enable_testing()
find_package(GTest) #find the gtest library  # will fail if not found

if(NOT GTest_FOUND)
  message(WARNING "GTest not found. Skipping targets.")
else()
  message(STATUS "GTest found.") # building tests 
  add_subdirectory(tests)
endif()

# adding google benchmark  
find_package(benchmark CONFIG)

if(NOT (benchmark_FOUND AND GTest_FOUND))
  message(WARNING "Google Benchmark not found. Skipping benchmark targets")
else()
  message(STATUS "Google Benchmark found.") #  building benchmarks
  add_subdirectory(benchmarks)
endif()

# examples subdirectory
add_subdirectory(examples)
