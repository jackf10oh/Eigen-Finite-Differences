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

#include<chrono>

namespace fdm{
  namespace solvers{ 

// Write Policies. i.e. write previous solution to cout, write to std::vector, write to CSV 
struct EmptySaver
{
  // no member data 
  EmptySaver()=default; 
  void saveSolution(Eigen::VectorXd sol){}; 
  void saveLastSolution(Eigen::VectorXd sol){}; 
};

struct LastSaver
{
  // no member data 
  LastSaver()=default; 
  void saveSolution(const Eigen::VectorXd& sol={}){}; 
  auto saveLastSolution(Eigen::VectorXd sol){ return sol; }; 
};

struct PrintSaver
{
  bool first_entry=true; 
  PrintSaver()=default; 
  void saveSolution(const Eigen::VectorXd& sol)
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
  void saveLastSolution(const Eigen::VectorXd& sol)
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

template< typename Units = std::chrono::milliseconds>
struct TimerSaver
{
  // single member data :-) 
  std::chrono::time_point<std::chrono::system_clock> time_started = std::chrono::system_clock::now(); 

  TimerSaver()=default; 

  void saveSolution(const Eigen::VectorXd& sol={}){}; 

  auto saveLastSolution(const Eigen::VectorXd& sol={})
  { 
    return std::chrono::duration_cast<Units>(std::chrono::system_clock::now() - time_started); 
  }; 
}; 

  } // end namespace solvers
} // end namespace fdm  

#endif // WritePolicies.hpp 