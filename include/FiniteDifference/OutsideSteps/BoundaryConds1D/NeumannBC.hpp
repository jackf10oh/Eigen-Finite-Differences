// NeumannBC.hpp 
//
// neumann boundary condtions of the form 
// Ux = c for some value c
//
// JAF 12/8/2025

#ifndef NEUMANNBCS_H
#define NEUMANNBCS_H 

#include "../../Utilities/FornbergCalc.hpp"
#include "BCPair.hpp" 

namespace fdm{
  namespace osteps{

class NeumannBC
{
  public:  
    // member data 
    double boundary_flux;

  public:
    // Constructors + destructor ---------------------------------------------
    NeumannBC(double val_init=0.0) : boundary_flux(val_init){}; 
    NeumannBC(const NeumannBC& other) : boundary_flux(other.boundary_flux){}; 
    // destructor
    ~NeumannBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    void SetStencilL(double t, const SharedConstMesh1D& mesh, fdm::Matrix& Mat) const 
    {
      Mat.topRows(1) *= 0;
      // first order derivative approximation 
      double h = mesh->operator[](1) - mesh->operator[](0);  
      Mat.coeffRef(0,0)= -1.0/h;
      Mat.coeffRef(0,1)=  1.0/h;
    }; 
    void SetStencilR(double t, const SharedConstMesh1D& mesh, fdm::Matrix& Mat) const 
    {
      Mat.bottomRows(1) *= 0; 
      // first order derivative approximation 
      double h = mesh->operator[](mesh->size()-1) - mesh->operator[](mesh->size()-2);  
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-2)= -1.0/h;
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=  1.0/h;
    };

    void SetImpSolL(double t, const SharedConstMesh1D& mesh, fdm::StridedRef Sol) const 
    {Sol[0] = boundary_flux;};
    void SetImpSolR(double t, const SharedConstMesh1D& mesh, fdm::StridedRef Sol) const 
    {Sol[Sol.size()-1] = boundary_flux;};
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(double t, const SharedConstMesh1D& mesh, fdm::StridedRef Sol) const 
    { 
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights((*mesh)[0], mesh->cbegin(), mesh->cbegin()+3, 1); 

      // solve the equation Flux = W[0]*S[0] + W[1]*S[1] + W[2]*S[2] 
      // for the target value S[0] 
      double target = boundary_flux; 
      target -= weights[1]*Sol[1]; 
      target -= weights[2]*Sol[2];
      target /= weights[0]; 

      // assign to Sol reference
      Sol[0] = target;  
      // void return type
    };
    void SetSolR(double t, const SharedConstMesh1D& mesh, fdm::StridedRef Sol) const 
    {
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights((*mesh)[mesh->size()-1], mesh->cend()-3, mesh->cend(), 1); 

      // solve the equation Flux = W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] 
      // for the target value S[N-1] 
      double target = boundary_flux; 
      target -= weights[0]*Sol[Sol.size()-3]; 
      target -= weights[1]*Sol[Sol.size()-2]; 
      target /= weights[2]; 

      // assign to Sol reference
      Sol[Sol.size()-1] = target;  
      // void return type
    };
};

  } // end namespace osteps
} // end namespace fdm 

#endif // NeumannBC.hpp