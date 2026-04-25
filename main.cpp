// main.cpp
//
//
//
// JAF 12/8/2025

#include<FiniteDifference/Diffops/Traits.hpp> 
#include<FiniteDifference/Mesh.hpp> 
#define EIGEN_SPARSEMATRIXBASE_PLUGIN <FiniteDifference/EigenFdmPlugin.hpp> 

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<memory>
#include<FiniteDifference/All.hpp> 
#include<FiniteDifference/Utilities/PrintVec.hpp> 
#include<Eigen/Dense>
#include<Eigen/Core>
#include<Eigen/SparseCore> // macro plugin takes effect. 
 
#include<FiniteDifference/Diffops/PartialDerivBase.hpp> 
#include<FiniteDifference/Diffops/NthPartialDeriv.hpp> 
#include<FiniteDifference/Diffops/NwiseUnaryOp.hpp> 
#include<FiniteDifference/Diffops/NwiseBinaryOp.hpp> 
#include<FiniteDifference/Diffops/EigenEvaluator.hpp> 

#include<FiniteDifference/Coeffs/CoeffBase.hpp>
#include<FiniteDifference/Coeffs/CoeffProduct.hpp> 
#include<FiniteDifference/Coeffs/AutonomousCoeff.hpp>
  
using namespace fdm; 

using std::endl, std::cout; 

int main()
{
  auto my_mesh = fdm::make_Mesh(2); 
  my_mesh->getAxis(0) = Eigen::VectorXd::LinSpaced(7,0.0,6.0); 
  my_mesh->getAxis(1) = Eigen::VectorXd::LinSpaced(7,0.0,6.0); 
  utils::print_vec(my_mesh->getAxis(0)); 

  linops::NthPartialDeriv<1,1> my_deriv; 

  linops::AutonomousCoeff my_coeff = [](double x, double y){ return x*x + 2*y; }; 

  auto mult = my_coeff * my_deriv; 
  mult.setMesh(my_mesh); 
  fdm::CSRMatrix result = mult; 
  cout << result << endl; 

  // auto expr = -my_deriv + 4.0 * my_deriv - my_deriv + 3 * my_deriv; 
  // expr.setMesh(my_mesh); 
  // fdm::CSRMatrix result = expr; 
  // cout << result << endl; 

};


// using MyArray = Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>; 
// MyArray weights(10); 
// weights << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;

// Eigen::Map<MyArray> map(weights.data(), 5); 
// cout << map << endl; 

// Eigen::Map<MyArray> map2(weights.data()+5, 5); 
// cout << map2 << endl; 

// using Mapped = Eigen::Map<Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>>;
// Mapped(weights.data(), 5) = Mapped(weights.data()+5, 5); 
// cout << weights << endl;  


// Future interfaces inside of diffops 

// Keep Mesh1D in the interface ? 
// adding a modulo + division to every InnerIterator might be costly ...

// template<int order, int dir, template<int minNodes> class NodeSelector>
// class NthDerivBase
// {
//   const Mesh* mesh_observed; 
//   std::weak_ptr<const Mesh> handle_observed; 
// }; 

// template<int order>
// class NodeSelector
// {
//   // constructor 
//   template<int dir>
//   NodeSelector(const NthDerivBase<order, dir, NodeSelector>& diffop)
//   {

//   }

//   // Fornberg algo needs order+1 nodes for finite difference weights 
//   static constexpr int numNodesMin = order+1; 
//   // possibly more nodes...  
//   static constexpr int numNodesMax = numNodesMin; 
//   // runtime int that lands in [numNodesMin, numNodesMax]
//   int numNodes; 

//   // offset where to write the innerIdxs into innerIndexPtr(), fornberg weights into valuePtr()
//   int ptr_offset;  
  
//   // Array of idxs that will be written into innerIndexPtr()[ptr_offset] 
//   std::array<std::size_t, numNodesMax> innerIdxs; 

//   // Array of doubles that will given to fornberg algorithm  
//   std::array<double, numNodesMax> nodes; 
// }; 

// // evaluator forward declaration ---------------- 
// template<int order, int dir, template<int minNodes> class NodeSelector, bool isNested>
// class evaluator{}; 

// // evaluator when it is the root of an expression 
// template<int order, int dir, template<int minNodes> class NodeSelector>
// struct evaluator<order,dir,NodeSelector, true>
// {
//   using Calculator = fdm::utils::FornArrayCalc<NodeSelector<order>::numNodesMax, order>; 
//   // Constructor  
//   evaluator(const NthDerivBase<order,dir,NodeSelector>& xpr)

//   // Member Functions ------------------------------------- 
//   void reseatMap(const Calculator& c, const int& numNodes)
//   {
//     new (&m_vals) = Eigen::::Map< Eigen::VectorXd>( c.getArray().data() + order * numNodes, numNodes);     
//   }

//   // Member Data --------------------------------------------  
//   // holds nodes + innerIdxs 
//   NodeSelector<order> m_selector; 

//   // fornberg calculator. use memory on the stack. 
//   Calculator m_calc; 

//   // Eigen expression of m_calc's weights 
//   Eigen::::Map< Eigen::VectorXd> m_vals; 

// }; 

// // evaluator when its nested in a parent evaluator 
// template<int order, int dir, template<int minNodes> class NodeSelector>
// struct evaluator<order,dir,NodeSelector, false>
// {
//   // Constructor: Use Parent's Fornberg Calc  
//   evaluator( 
//     const NthDerivBase<order,dir,NodeSelector>& xpr
//     const fdm::utils::FornArrayCalc<NodeSelector<order>::numNodesMax, order>& parent_calc, 
//     const int& numNodes)
//   {
//     // initialize reference to NthDerivBase<...> before this call? 
//     reseatMap(parent_calc, numNodes)
//   }

//   // Member Data ------------------------ 
//   // doesn't need it's own NodeSelector or FornArrayCalc

//   // Member Functions ------------------------------------- 
//   void reseatMap(const Calculator& c, const int& numNodes)
//   {
//     new (&m_vals) = Eigen::::Map< Eigen::VectorXd>( c.getArray().data() + order * numNodes, numNodes);     
//   }
  
//   // TODO Eigen expression of parent_calc's weights 
//   Eigen::Map<Eigen::VectorXd> m_vals; 
// }; 

// // --------------------------------------------

