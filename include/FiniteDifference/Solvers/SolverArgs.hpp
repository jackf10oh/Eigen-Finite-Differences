// SolverArgs.hpp
//
// P.O.D. class containing mesh of time, mesh of space, and initial conditions. 
// Initial conditions are specified by first N solution at first N entries of time. 
//
// JAF 3/4/2026 

#ifndef FDM_SOLVERS_SOLVERARGS_H
#define FDM_SOLVERS_SOLVERARGS_H

#include<memory>
#include<vector>
#include<Eigen/Core> 
#include "../Types.hpp"
#include "../Mesh.hpp"


namespace fdm{
  namespace solvers{ 

template<class M, class Container>
struct SolverArgs
{
  // Mesh1D or MeshXD the PDE operates on 
  std::shared_ptr<M> mesh; 

  // list of times the solver marches through 
  std::shared_ptr<const Container> times; 
  
  // first N solution values. 
  // ! has to be atleast >= maxOrder + 1. 
  // where maxOrder is the highest order in the LHS time derivatives expression (texprs)  
  // defaulted to empty so it can be assigned later... 
  std::vector<fdm::Vector> initialConditions = {};
};

// CTAD guideline ... 
template<class M, class C>
SolverArgs(std::shared_ptr<M>, std::shared_ptr<const C>, std::vector<fdm::Vector>)
  ->SolverArgs<M, C>; 

template<class M, class C>
SolverArgs(std::shared_ptr<M>, std::shared_ptr<const C>)
  ->SolverArgs<M, C>; 

  } // end namespace solvers
} // end namespace fdm 

#endif /// SolverArgs.hpp 