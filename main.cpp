// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include<LinOps/All.hpp>  
// #include<OutsideSteps/All.hpp> 
// #include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 
#include<TExprs/All.hpp> 
// #include<Solvers/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/FornbergArrayCalc.hpp> 

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  constexpr std::size_t order = 2; 
  constexpr std::size_t num_nodes = 3; 
  double x_bar = 1.0; 

  FornArrayCalc<num_nodes, order> calc;
  FornCalc calc2(num_nodes, order);
  
  LinOps::Mesh1D_SPtr_t m = LinOps::make_mesh(0.0, num_nodes - 1, num_nodes); 

  calc.calculate(x_bar, m->cbegin(), m->cend()); 
  calc2.Calculate(x_bar, m->cbegin(), m->cend(), order); 
 
  // print_vec(calc.getArray(), "array calc"); 
  // print_vec(calc2.m_arr, "vector calc"); 

  TExprs::NthTimeDeriv<2> Utt{};
  auto c_func =  [](double x){ return x*x; }; 
  auto c = LinOps::AutonomousCoeff(c_func); 

  auto mult = 2.0 * Utt; 
  auto mult02 = c * Utt; 

  auto lam = [&](auto&& arg){ cout << "arg is lval? " << std::is_lvalue_reference<decltype(arg)>::value << "\n"; }; 

  lam(Utt); 
  lam(mult); 
  lam(mult02); 

  auto lam02 = [&](auto&& arg){ cout << "store by lval? " << std::is_lvalue_reference<typename TExprs::traits::Storage<decltype(arg)>::type>::value << "\n"; };
  lam02(Utt); 
  lam02(mult); 
  lam02(mult02); 

  cout << "-----------------------\n"; 
  cout << "1." << endl; 
  auto lam03 = [&](auto&& arg){ cout << "numNodes: " << arg.numNodes << " size of times: " << arg.getStoredTimes().size() << " size of sols: " << arg.getStoredSolutions().size() << "\n"; }; 
  auto x1 = TExprs::make_Executor<4>(Utt); 
  lam03(x1); 

  cout << "2." << endl; 
  TExprs::Executor<decltype(mult),6> x2(mult); 
  lam03(x2); 

  cout << "3." << endl; 
  auto x3 = TExprs::make_Executor<5>(mult02); 
  lam03(x3); 


  using TExprs::NthTimeDeriv; 
  // auto messy = NthTimeDeriv<7>{} - NthTimeDeriv<4>{} + NthTimeDeriv<2>{}; 
  auto messy = NthTimeDeriv<7>{} - NthTimeDeriv<4>{} + NthTimeDeriv<2>{}; 

  auto tup = messy.toTuple(); 
  cout << "tuple[0] is lval? " << std::is_lvalue_reference<std::tuple_element_t<0,decltype(tup)>>::value << "\n"; 
  cout << "tuple[1] is lval? " << std::is_lvalue_reference<std::tuple_element_t<1,decltype(tup)>>::value << "\n"; 
  cout << "tuple[2] is lval? " << std::is_lvalue_reference<std::tuple_element_t<2,decltype(tup)>>::value << "\n"; 
  
  cout << "finish -----------------------\n"; 

  cout << "-----------------------\n"; 

  auto result01 = TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_double>(mult.toTuple()); 
  auto result02 = TExprs::traits::filter_tup<TExprs::traits::coeffat_returns_other>(mult.toTuple()); 
  cout << "size of tuple 01: " << std::tuple_size_v<decltype(result01)> << "\n"; 
  cout << "size of tuple 02: " << std::tuple_size_v<decltype(result02)> << "\n"; 

  cout << "tuple[0] is lval? " << std::is_lvalue_reference<std::tuple_element<0,decltype(mult.toTuple())>::type>::value << "\n"; 
  cout << "tuple[0] is lval? " << std::is_lvalue_reference<std::tuple_element<0,decltype(result01)>::type>::value << "\n"; 

  auto x = TExprs::make_Executor(messy); 
  cout << x.numNodes << endl; 
};