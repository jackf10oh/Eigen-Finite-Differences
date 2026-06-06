// NeumannBC.hpp 
//
// neumann boundary condtions of the form 
// Ux = c for some value c
//
// JAF 12/8/2025

#ifndef FORNFDM_OSTEPS_NEUMANNBC_H
#define FORNFDM_OSTEPS_NEUMANNBC_H

#include "../../Utilities/FornbergStackCalc.hpp"
#include "BCPair.hpp" 

namespace fornfdm{
  namespace osteps{

class NeumannBC
{
  public:  
    // member data 
    fornfdm::Scalar flux;

  public:
    // Constructors + destructor ---------------------------------------------
    NeumannBC(fornfdm::Scalar val_init=0.0) : flux(val_init){}; 
    NeumannBC(const NeumannBC& other) : flux(other.flux){}; 
    // destructor
    ~NeumannBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fornfdm stencil matrix
    BoundaryRow getTopRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    {
      // first order derivative approximation. TODO use 3 node stencil here 
      fornfdm::Scalar h = axis[1] - axis[0];
      return BoundaryRow{{-1.0/h, 1.0/h, 0}, 2};  
    }; 
    BoundaryRow getBottomRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    {
      // first order derivative approximation 
      fornfdm::Real h = axis[axis.size()-1] - axis[axis.size()-2];  
      return BoundaryRow{{-1.0/h, 1.0/h, 0}, 2};  
    };

    void SetImpSolL(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {Sol[0] = flux;};
    void SetImpSolR(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {Sol[Sol.size()-1] = flux;};
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    { 
      assert((Sol.size()>=3 || mesh.size()>=3) && "Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fornfdm::utils::FornbergStackCalc<3,1> calc; 

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      calc.calculate(mesh[0], mesh.cbegin(), mesh.cbegin()+3); 

      // solve the equation Flux = W[0]*S[0] + W[1]*S[1] + W[2]*S[2] 
      // for the target value S[0] 
      fornfdm::Scalar target = flux; 
      target -= calc.getArray()[4]*Sol[1]; 
      target -= calc.getArray()[5]*Sol[2];
      target /= calc.getArray()[3]; 

      // assign to Sol reference
      Sol[0] = target;  
    };
    void SetSolR(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef Sol) const 
    {
      assert((Sol.size()>=3 || mesh.size()>=3) && "Vector or Axis size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fornfdm::utils::FornbergStackCalc<3,1> calc; 

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      calc.calculate(mesh[mesh.size()-1], mesh.cend()-3, mesh.cend()); 

      // solve the equation Flux = W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] 
      // for the target value S[N-1] 
      fornfdm::Scalar target = flux; 
      target -= calc.getArray()[3]*Sol[Sol.size()-3]; 
      target -= calc.getArray()[4]*Sol[Sol.size()-2]; 
      target /= calc.getArray()[5]; 

      // assign to Sol reference
      Sol[Sol.size()-1] = target;  
    };
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // NeumannBC.hpp