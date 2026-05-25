// DirichletBC.hpp 
//
// dirichlet boundary value conditions of the form 
// U(0) = c for some constant c 
//
// JAF 12/8/2025

#ifndef FORNFDM_OSTEPS_DIRICHLETBC_H
#define FORNFDM_OSTEPS_DIRICHLETBC_H

#include "BCPair.hpp" 

namespace fornfdm{
  namespace osteps{

class DirichletBC
{
  public:  
    // member data 
    fornfdm::Scalar value;

  public:
    // Constructors ---------------------------------------------
    DirichletBC(fornfdm::Scalar val_init=0.0) : value(val_init){}; 
    DirichletBC(const DirichletBC& other)=default; 
    ~DirichletBC()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fornfdm stencil matrix
    void SetStencilL(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::CSRMatrix& Mat) const 
    {
      Mat.topRows(1) *= 0; Mat.coeffRef(0,0)=1;
    }
    void SetStencilR(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::CSRMatrix& Mat) const 
    {
      Mat.bottomRows(1) *= 0; Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=1;
    }

    // change the first/last (left/right boundary) entry of a vector to implicit solution   
    void SetImpSolL(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {Sol[0] = value;}
    void SetImpSolR(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {Sol[Sol.size()-1] = value;}
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    { Sol[0] = value;}
    void SetSolR(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {Sol[Sol.size()-1] = value;}
};


  } // end namespace osteps
} // end namespace fornfdm 

#endif // DirichletBC.hpp