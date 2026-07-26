cmake_minimum_required(VERSION 3.16)
project(fornfdm LANGUAGES CXX)
add_library(fornfdm INTERFACE)

# configure ctest to find test added by subdirectories
include(FetchContent)

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
option(USE_OPENMP "Configure Cmake to build GTest suite" OFF) 
if(BUILD_TESTING OR BUILD_BENCHMARKS)
  enable_testing()
  # find_package(GTest) #find the gtest library  # will fail if not found
  FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG v1.14.0
  )
  FetchContent_MakeAvailable(googletest)
  if(TARGET GTest::gtest)
    message(STATUS "GTest found.") # building tests 
    add_subdirectory(tests)
  else()
    message(WARNING "GTest not found.")
  endif()
endif()

# adding google benchmark  
if(BUILD_BENCHMARKS)
  find_package(benchmark CONFIG)
  if(benchmark_FOUND)
    message(STATUS "Google Benchmark found.") #  building benchmarks
    add_subdirectory(benchmarks)
  else()
    message(WARNING "Google Benchmark not found. Skipping benchmark targets")
  endif()
endif()

# examples subdirectory
add_subdirectory(examples)

# cmake .. `
#   -G "MinGW Makefiles" `
#   -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe `
#   -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe `
#   -DEigen3_DIR=C:/eigen-install/share/eigen3/cmake `