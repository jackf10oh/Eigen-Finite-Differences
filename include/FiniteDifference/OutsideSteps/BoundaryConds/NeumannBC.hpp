// NeumannBC.hpp 
//
// neumann boundary condtions of the form 
// Ux = c for some value c
//
// JAF 12/8/2025

#ifndef FDM_OSTEPS_NEUMANNBC_H
#define FDM_OSTEPS_NEUMANNBC_H

#include "../../Utilities/FornbergStackCalc.hpp"
#include "BCPair.hpp" 

namespace fdm{
  namespace osteps{

class NeumannBC
{
  public:  
    // member data 
    fdm::Scalar boundary_flux;

  public:
    // Constructors + destructor ---------------------------------------------
    NeumannBC(fdm::Scalar val_init=0.0) : boundary_flux(val_init){}; 
    NeumannBC(const NeumannBC& other) : boundary_flux(other.boundary_flux){}; 
    // destructor
    ~NeumannBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    void SetStencilL(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.topRows(1) *= 0;
      // first order derivative approximation 
      fdm::Scalar h = mesh[1] - mesh[0];  
      Mat.coeffRef(0,0)= -1.0/h;
      Mat.coeffRef(0,1)=  1.0/h;
    }; 
    void SetStencilR(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.bottomRows(1) *= 0; 
      // first order derivative approximation 
      fdm::Real h = mesh[mesh.size()-1] - mesh[mesh.size()-2];  
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-2)= -1.0/h;
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=  1.0/h;
    };

    void SetImpSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[0] = boundary_flux;};
    void SetImpSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[Sol.size()-1] = boundary_flux;};
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    { 
      assert((Sol.size()<3 || mesh.size()<3) && "Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornbergStackCalc<3,1> calc; 

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      calc.calculate(mesh[0], mesh.cbegin(), mesh.cbegin()+3); 

      // solve the equation Flux = W[0]*S[0] + W[1]*S[1] + W[2]*S[2] 
      // for the target value S[0] 
      fdm::Scalar target = boundary_flux; 
      target -= calc.getArray()[4]*Sol[1]; 
      target -= calc.getArray()[5]*Sol[2];
      target /= calc.getArray()[3]; 

      // assign to Sol reference
      Sol[0] = target;  
    };
    void SetSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {
      assert((Sol.size()<3 || mesh.size()<3) && "Vector or Axis size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornbergStackCalc<3,1> calc; 

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      calc.calculate(mesh[mesh.size()-1], mesh.cend()-3, mesh.cend()); 

      // solve the equation Flux = W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] 
      // for the target value S[N-1] 
      fdm::Scalar target = boundary_flux; 
      target -= calc.getArray()[3]*Sol[Sol.size()-3]; 
      target -= calc.getArray()[4]*Sol[Sol.size()-2]; 
      target /= calc.getArray()[5]; 

      // assign to Sol reference
      Sol[Sol.size()-1] = target;  
    };
};

  } // end namespace osteps
} // end namespace fdm 

#endif // NeumannBC.hpp