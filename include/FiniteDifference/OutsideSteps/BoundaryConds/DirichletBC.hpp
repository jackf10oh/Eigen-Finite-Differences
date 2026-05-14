// DirichletBC.hpp 
//
// dirichlet boundary value conditions of the form 
// U(0) = c for some constant c 
//
// JAF 12/8/2025

#ifndef FDM_OSTEPS_DIRICHLETBC_H
#define FDM_OSTEPS_DIRICHLETBC_H

#include "BCPair.hpp" 

namespace fdm{
  namespace osteps{

class DirichletBC
{
  public:  
    // member data 
    fdm::Scalar boundary_val;

  public:
    // Constructors ---------------------------------------------
    DirichletBC(fdm::Scalar val_init=0.0) : boundary_val(val_init){}; 
    DirichletBC(const DirichletBC& other)=default; 
    ~DirichletBC()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    void SetStencilL(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.topRows(1) *= 0; Mat.coeffRef(0,0)=1;
    }
    void SetStencilR(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.bottomRows(1) *= 0; Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=1;
    }

    // change the first/last (left/right boundary) entry of a vector to implicit solution   
    void SetImpSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[0] = boundary_val;}
    void SetImpSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[Sol.size()-1] = boundary_val;}
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    { Sol[0] = boundary_val;}
    void SetSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[Sol.size()-1] = boundary_val;}
};


  } // end namespace osteps
} // end namespace fdm 

#endif // DirichletBC.hpp