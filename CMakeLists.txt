cmake_minimum_required(VERSION 3.16)
project(FornFdm LANGUAGES CXX)
add_library(FornFdm INTERFACE)

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
    target_link_libraries(FornFdm INTERFACE OpenMP::OpenMP_CXX)
  endif()
else()
  message(STATUS "Cmake args disabling OpenMP...")
endif()

# adding Eigenlib to compiler
find_package(Eigen3 REQUIRED)

# header only library ===============================
target_include_directories(FornFdm 
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(FornFdm
  INTERFACE Eigen3::Eigen
)

# adding google test
enable_testing()
find_package(GTest) #find the gtest library  # will fail if not found

# adding google benchmark  
# find_package(benchmark) # Google Benchmark will fail if not found
# find_package(benchmark REQUIRED) # Google Benchmark will fail if not found

# if(NOT (benchmark_FOUND AND GTest_FOUND))
#   message(WARNING "Google Benchmark not found. Skipping benchmark targets")
# else()
#   message(STATUS "Google Benchmark found.") #  building benchmarks
#   # add_executable(benchmarks
#   #   benchmarks/benchmarks.cpp
#   # ) 
#   # target_link_libraries(benchmarks
#   #   PRIVATE LinOps 
#   #   PRIVATE benchmark::benchmark
#   #   PRIVATE GTest::gtest 
#   #   PRIVATE GTest::gtest_main 
#   #   PRIVATE pthread # required by GTest
#   # )

#   # if(USE_OPENMP)
#   #   target_link_libraries(benchmarks
#   #     PRIVATE OpenMP::OpenMP_CXX
#   #   )
#   # endif()

# endif()

if(NOT GTest_FOUND)
  message(WARNING "GTest not found. Skipping targets.")
else()
  message(STATUS "GTest found.") # building tests 
  add_subdirectory(tests)
endif()

add_executable(main 
  main.cpp
)

target_link_libraries(main
  PRIVATE FornFdm
)

# on windows with mingw 
# cmake .. `
#   -G "MinGW Makefiles" `
#   -DCMAKE_C_COMPILER=C:/msys64/mingw64/bin/gcc.exe `
#   -DCMAKE_CXX_COMPILER=C:/msys64/mingw64/bin/g++.exe `
#   -DEigen3_DIR=C:/eigen-install/share/eigen3/cmake `
#   -DGTest_DIR=C:\gtest-mingw\lib\cmake\GTest `
#   -Dpybind11_DIR=C:/Users/jackf/anaconda3/Lib/site-packages/pybind11/share/cmake/pybind11
