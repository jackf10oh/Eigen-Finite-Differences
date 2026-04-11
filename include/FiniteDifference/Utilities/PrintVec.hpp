// PrintVec.hpp 
//
// prints vectors in a pythonic way: 
// i.e. comment: [x0, x1, ..., xn] 
//
// JAF 1/4/2026 

#ifndef PRINTVEC_H
#define PRINTVEC_H

#include<iostream>
#include<Eigen/Core> 

// from container with .begin() .end() 
template<typename Cont>
void print_vec(const Cont& v, std::string comment="", bool new_line=true){
  // print comment 
  if(!comment.empty()) std::cout << comment << ": ";
  // print [ x0, x1, ..., xn] 
  std::cout << "[" << *std::for_each_n(v.cbegin(), v.size()-1, [](const auto& x){std::cout << x << ", ";}) << "]"; 
  // print new line 
  if(new_line) std::cout << "\n";   
};

// from iterators. 
template<typename Iter>
void print_vec(const Iter& start, const Iter& stop, std::string comment="", bool new_line=true){
  // print comment 
  if(!comment.empty()) std::cout << comment << ": ";
  // print [ x0, x1, ..., xn] 
  std::cout << "[" << *std::for_each_n(start, std::distance(start,stop)-1, [](const auto& x){std::cout << x << ", ";}) << "]"; 
  // print new line 
  if(new_line) std::cout << "\n";   
};

// for Eigen matrix 
template<typename Derived>
void print_mat(const Eigen::DenseBase<Derived>& A, std::string comment = "")
{
  // print comment 
  if(!comment.empty()) std::cout << "--------------" << comment << "--------------" << "\n";
  // opening brace 
  std::cout << "["; 
  // print rows up to last row with , << \n 
  for(std::size_t i=0; i<A.rows()-1; i++){
    print_vec(A.row(i), "", false); 
    std::cout << ",\n"; 
  }
  // print last row
  print_vec(A.row(A.rows()-1), "", false); 
  // closing brace 
  std::cout << "]\n\n";
}; 


// for an STL like container 
template<typename Cont2D, 
  typename = std::enable_if_t<
    !std::is_base_of<Eigen::DenseBase<Cont2D>, Cont2D>::value
  >
>
void print_mat(const Cont2D& A, std::string comment="")
{
  // print comment 
  if(!comment.empty()) std::cout << "--------------" << comment << "--------------" << "\n";
  // opening brace 
  std::cout << "["; 
  // print rows up to last row with , << \n 
  auto it_last_row = std::for_each_n(A.begin(), A.size()-1, 
      [](const auto& row){
          print_vec(row,"",false); std::cout<<","<< "\n";
      }
  ); 
  // print last row
  print_vec(*it_last_row, "", false); 
  // closing brace 
  std::cout << "]\n\n";
}

#endif // PrintVec.hpp