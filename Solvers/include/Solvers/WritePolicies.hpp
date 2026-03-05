// WritePolicies.hpp 
//
// Instructs a SolverClass how to handle Solutions U(0), U(1), ..., U(n-1), U(n) 
// after they have been calculated. 
// 
// only 2 points of variation. 
// 1. Handling intermediate solutions U(0), U(1), ..., U(n-1)
// 2. Handling the last solution U(n) 
//
// JAF 3/4/2025 

#ifndef WRITEPOLICY_H
#define WRITEPOLICY_H 

namespace Solvers{ 

// Write Policies. i.e. write previous solution to cout, write to std::vector, write to CSV 
struct EmptyWrite
{
  // no member data 
  EmptyWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol){}; 
  void ConsumeLastSolution(Eigen::VectorXd&& sol){}; 
};

struct FinalWrite
{
  // no member data 
  FinalWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol){}; 
  auto ConsumeLastSolution(Eigen::VectorXd&& sol){ return sol; }; 
};

struct PrintWrite
{
  bool first_entry=true; 
  PrintWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol)
  {
    if(first_entry){
      std::cout << "[";
      first_entry=false;
    }
    std::cout << "["; 
    auto it=sol.cbegin(); 
    auto end = std::prev(sol.cend()); 
    for(; it!=end; it++){
      std::cout << *it << ", ";
    }
    std::cout << *it << "],\n";  
  }; 
  void ConsumeLastSolution(Eigen::VectorXd&& sol)
  {
    std::cout << "["; 
    auto it=sol.cbegin(); 
    auto end = std::prev(sol.cend()); 
    for(; it!=end; it++){
      std::cout << *it << ", ";
    }
    std::cout << *it << "]]\n";  
  }; 
};

} // end namespace Solvers 

#endif // WritePolicies.hpp 