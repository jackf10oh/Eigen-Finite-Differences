// main.cpp
//
// JAF 12/8/2025

// workaround to allow expressions of different rows/cols to be added 
// before setMesh() sets their rows/cols to be equal
#define eigen_assert(x)
#include<FornFdm/All.hpp>
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <FornFdm/EigenFdmPlugin.hpp> 

#include<iostream>
#include<iomanip>
#include<FornFdm/Utilities/PrintVec.hpp> 
#include<FornFdm/Utilities/BumpFunc.hpp>
#include<Eigen/SparseCore> // macro plugin takes effect. 

using std::endl, std::cout; 

int main()
{
  using namespace fornfdm; 

  // meshes are passed to fornfdm classes by shared_ptr

  // 1 dimensionsal mesh 
  std::shared_ptr<fornfdm::Mesh> mesh_01 = fornfdm::make_Mesh(1); 

  // 2 dimensional mesh
  auto mesh_02 = fornfdm::make_Mesh(2);

  // 2 axes with copies of the same value on construction
  auto mesh_03 = fornfdm::make_Mesh(fornfdm::linspaced(11,0.0,10.0), 2); 

  // 3 axes each with different initials 
  auto mesh_04 = make_Mesh(linspaced(11,0.0,10.0), linspaced(21, 0.0, 20.0), linspaced(101,-10.0,10.0)); 

  // getters 
  std::cout << "# of dims: " << mesh_04->numDims() << std::endl; 
  std::cout << "axis 0: " << mesh_04->getAxis(0).transpose() << std::endl;
  std::cout << "axis 0 size: " << mesh_04->sizeOfDim(0) << std::endl;
  std::cout << "axis 1: " << mesh_04->getAxis(1).transpose() << std::endl;
  std::cout << "axis 1 size: " << mesh_04->sizeOfDim(1) << std::endl;
  std::cout << "axis 1: " << mesh_04->getAxis(1).transpose() << std::endl;
  std::cout << "axis 1 size: " << mesh_04->sizeOfDim(1) << std::endl;

  // utilities 
  std::cout << "full sizes product: " << mesh_04->sizesProduct() << std::endl; 
  std::cout << "middle product [0,2]: " << mesh_04->sizesMiddleProduct(0,2) << std::endl; 

  // swappable. (asserts # of dimensions are equal)
  mesh_02->swap(*mesh_03); 


  std::cout << "\n New start ------------------------- \n" << std::endl; 
  fornfdm::Coordinate c{0.0,1.0,2.0};
  cout << "coord: " << c << endl; 
};
