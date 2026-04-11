// SolverArgs.hpp
//
// P.O.D. class containing mesh of time, mesh of space, and initial conditions. 
// Initial conditions are specified by first N solution at first N entries of time. 
//
// JAF 3/4/2026 

#ifndef SOLVERARGS_H
#define SOLVERARGS_H

#include<LinOps/Mesh.hpp> 

namespace Solvers{ 

template<typename ANYMESH_SHAREDPTR_T = LinOps::SharedConstMesh1D>
struct SolverArgs
{
  // shared_ptr to const Mesh1D or const MeshXD 
  ANYMESH_SHAREDPTR_T domain_mesh_ptr; 

  // shared_ptr to const Mesh1D 
  LinOps::SharedConstMesh1D time_mesh_ptr; 
  
  // List of underlying Eigen::VectorXd from Discretization1D (XD)  
  std::vector<Eigen::VectorXd> ICs;
}; 

// CTAD guideline ... 
template<typename M>
SolverArgs(M,LinOps::SharedConstMesh1D, std::vector<Eigen::VectorXd>)
  ->SolverArgs<M>; 

} // end namespace Solvers 

#endif /// SolverArgs.hpp 