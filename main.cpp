// main.cpp
//
//
//
// JAF 12/8/2025

#include<FiniteDifference/Traits.hpp> 
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
#include<Eigen/SparseCore> // macro plugin takes effect. 

#include<FiniteDifference/Diffops/NodeSelector.hpp>
  
using namespace fdm; 

using std::endl, std::cout; 

int main()
{
  auto my_mesh = make_Mesh1D(0.0,10.0,11); 

  linops::NthDerivOp<1> my_deriv;
  my_deriv.setMesh(my_mesh); 

  fdm::CSRMatrix A = my_deriv.asMatrix(); 
  cout << A << endl; 

  cout << "nnz A: " << A.nonZeros() << " estimate: " << linops::internal::NodeSelector<2>::sumNodesPerRow(*my_mesh) << endl; 

  for(auto i=0; i<A.rows(); i++){
    auto start = A.outerIndexPtr()[i];
    auto end = A.outerIndexPtr()[i+1];
    cout << "row " << i <<": " << (end-start) << "idxs: ";  
    for(auto j=start; j<end; ++j){
      cout << A.innerIndexPtr()[j] << ", "; 
    }
    cout << "\n";
    
    linops::internal::NodeSelector<2> test_selector(*my_mesh, i); 
    cout << test_selector.numNodesUsed << " idxs: "; 
    for(auto j=0; j<test_selector.numNodesUsed; ++j){
      cout << test_selector.nodeIndices[j] << ", "; 
    }
    cout << "\n"; 
  }

  bool success = true; 
  if(A.nonZeros() != linops::internal::NodeSelector<2>::sumNodesPerRow(*my_mesh)) success = false;

  for(auto i=0; i<A.rows(); i++){
    linops::internal::NodeSelector<2> test_selector(*my_mesh, i); 
    auto start = A.outerIndexPtr()[i];
    auto end = A.outerIndexPtr()[i+1];
    for(auto j=start; j<end; ++j){
      if(A.innerIndexPtr()[j] != test_selector.nodeIndices[j-start]) success = false;  
    }
    if((end - start) != test_selector.numNodesUsed) success = false; 
  }

  cout << "everything looks good? " << success << endl; 


};

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

