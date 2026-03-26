// DirichletBC.hpp 
//
// dirichlet boundary value conditions of the form 
// U(0) = c for some constant c 
//
// JAF 12/8/2025

#ifndef DIRICHLETBCS_H
#define DIRICHLETBCS_H 

#include "BCPair.hpp" 

namespace OSteps{

class DirichletBC
{
  public:  
    // member data 
    double boundary_val;

  public:
    // Constructors ---------------------------------------------
    DirichletBC(double val_init=0.0) : boundary_val(val_init){}; 
    DirichletBC(const DirichletBC& other)=default; 
    ~DirichletBC()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    void SetStencilL(double t, const LinOps::Mesh1D_SPtr_t& mesh, MatrixStorage_t& Mat) const 
    {
      Mat.topRows(1) *= 0; Mat.coeffRef(0,0)=1;
    }
    void SetStencilR(double t, const LinOps::Mesh1D_SPtr_t& mesh, MatrixStorage_t& Mat) const 
    {
      Mat.bottomRows(1) *= 0; Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=1;
    }

    // change the first/last (left/right boundary) entry of a vector to implicit solution   
    void SetImpSolL(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const 
    {Sol[0] = boundary_val;}
    void SetImpSolR(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const 
    {Sol[Sol.size()-1] = boundary_val;}
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const 
    { Sol[0] = boundary_val;}
    void SetSolR(double t, const LinOps::Mesh1D_SPtr_t& mesh, StridedRef_t Sol) const 
    {Sol[Sol.size()-1] = boundary_val;}
};

} // end namespace OSteps 

#endif // DirichletBC.hpp