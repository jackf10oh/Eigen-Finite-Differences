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

using namespace fornfdm; 

using std::endl, std::cout; 

int main()
{
  std::cout << std::setprecision(4);
  auto mesh = fornfdm::make_Mesh(fornfdm::linspaced(5,0.0,4.0),2); 

  fornfdm::CSRMatrix mat; 
  mat.resize(25, 25); 
  mat = Eigen::MatrixXd::Random(25,25).sparseView();
  fornfdm::Vector sol = fornfdm::make_Discretization(mesh, 0.0); 
  // cout << mat.toDense() << endl;

  auto bc1 = osteps::DirichletBC(1.0); 
  auto bc2 = osteps::DirichletBC(2.0);
  // auto bc3 = osteps::RobinBC(5.0,1.0,30.0); 

  using osteps::BCPair;
  osteps::BCList bcs(BCPair(bc1,bc1), BCPair(bc2,bc2));

  auto t = osteps::make_time(); 
  auto ctx = osteps::make_context(mesh); 
  // bcs.MatBeforeStep<osteps::StepType::Implicit>(mat,t,ctx);
  // cout << mat.toDense() << endl;

  bcs.VecAfterStep<osteps::StepType::Explicit>(sol, t, ctx);
  cout << sol << endl;
};

// for(block in big blocks)
// first blocks is set by lower blocks but masked by this row iterators
// recursively goes into lower dimension

// void setStencilNoFill(fornfdm::CSRMatrix& mat, std::size_t offset)
// {
//   if(ith_dimension == 0)
//   {
//     // use a row iterator to fill the top rows. 
//     // use a row iterator to fill the bottom rows.
//     // return/goto parent dimension
//   }
//   else
//   {
//     // setStencilFill() using THIS dimensions row iterator
//     // fill middle blocks with setStencilNoFill() and a certain offset 
//     // setStencilFill() using THIS dimensions row iterator 
//   }
// }

// void setStencilFill()
// {
//   if(ith_dimension == 0)
//   {
//     // use a row iterator to fill the top rows. 
//     // use a row iterator to fill middle rows. NEED a row iterator from parent dimension
//     // use a row iterator to fill the bottom rows. 
//     // How to determine parent dimension ?????
//   }
//   else
//   {
//     // setStencilFill() top block but THIS as row iterator
//     // setStencillFill() block but PARENT as row iterator
//     // setStencilFill() bottom block but THIS as row iterator
//   }
// }