// main.cpp
//
//
//
// JAF 12/8/2025

#include<Traits.hpp> 
#include<FiniteDifference/Mesh.hpp> 
#define EIGEN_SPARSEMATRIXBASE_PLUGIN "EigenBasePlugin.hpp" 

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<memory>
// #include<FiniteDifference/All.hpp> 
// #include<FiniteDifference/Mesh.hpp> 
// #include<Eigen/Dense>
// #include<FiniteDifference/Utilities/PrintVec.hpp> 
// #include<FiniteDifference/Utilities/BumpFunc.hpp> 

#include<Eigen/SparseCore> // macro plugin takes effect. 
#include<Foo.hpp> 
  
using namespace fdm; 

using std::endl, std::cout; 
int main()
{

  fdm::Matrix A, B; 
  A.resize(10,10); 
  B.resize(10,10); 


  // // testing the plugin 
  auto mesh = fdm::make_Mesh(); 
  Foo my_linop; 
  my_linop.resize(10,10); 

  auto messy = A + my_linop; 
  messy.setMesh(mesh); 
  messy.setTime(3.0); 

  cout << endl << "--------------------" << endl; 
  cout << "messy's mesh is null? " << (messy.getMesh()==nullptr) << endl; 
  cout << "A's mesh is null? " << (A.getMesh()==nullptr) << endl; 
  cout << "my_linops's mesh is null? " << (my_linop.getMesh()==nullptr) << endl; 
  cout << "rhs mesh is null? " << (messy.rhs().derived().getMesh()==nullptr) << endl; 

  cout << endl << "--------------------" << endl; 
  cout << "messy's time: " << messy.getTime() << endl; 
  cout << "A's time: " << A.getTime() << endl; 
  cout << "my_linop's time: " << my_linop.getTime() << endl; 
  cout << "rhs time: " << messy.rhs().getTime() << endl; 

  cout << endl << "--------------------" << endl; 
  cout << "messy is linop? " << traits::is_linop<decltype(messy)>::value << endl; 
  cout << "A is linop? " << traits::is_linop<decltype(A)>::value << endl; 
  cout << "my_linop is linop? " << traits::is_linop<decltype(my_linop)>::value << endl; 

  // cout << "A storage: " <<

  cout << endl << "--------------------" << endl; 
  cout << "manual func ptr == " << (&Foo::setTime == &Foo::Base::setTime) << endl; 
  cout << "decltype is time dep? " << traits::is_time_dep<decltype(my_linop)>::value << endl;
  cout << "pure type is_time_dep? " << traits::is_time_dep<Foo>::value << endl; 
  cout << "A is_time_dep? " << traits::is_time_dep<decltype(A)>::value << endl; 
  cout << "messy is_time_dep? " << traits::is_time_dep<decltype(messy)>::value << endl; 

  // cout << "pure type is linop? " << traits::is_linop<Foo>::value << endl; 


  // // const auto& left = messy.lhs(); 
  // auto& left = messy.lhs().const_cast_derived(); 

  // // left.setMesh(mesh); 

  // cout << endl << "--------------------" << endl; 
  // cout << "messy::RhsNested is lval? " << std::is_lvalue_reference<decltype(messy)::RhsNested>::value << endl; 
  // cout << "messy::RhsNested is const? " << std::is_const<decltype(messy)::RhsNested>::value << endl; 
};
