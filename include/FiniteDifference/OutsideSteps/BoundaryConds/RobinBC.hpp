// Robin.hpp 
//
// robin(third type) boundary conditions 
// alpha*U + beta*Ux = target_val
//
// JAF 12/8/2025

#ifndef FDM_OSTEPS_ROBINBC_H
#define FDM_OSTEPS_ROBINBC_H

#include "../../Utilities/FornbergCalc.hpp"
#include "BCPair.hpp"

namespace fdm{
  namespace osteps{

class RobinBC
{
  public:  
    // member data 
    fdm::Scalar boundary_target;
    fdm::Scalar val_coeff;
    fdm::Scalar deriv_coeff;

  public:
    // Constructors + Destructor ---------------------------------------------
    // a*U + b*Ux = target
    RobinBC(fdm::Scalar a=1.0, fdm::Scalar b=0.0, fdm::Scalar target=0.0) 
      : val_coeff(a), deriv_coeff(b), boundary_target(target)
    {}; 
    // copy 
    RobinBC(const RobinBC& other)=default; 
    // destructors
    virtual ~RobinBC()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    void SetStencilL(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.topRows(1) *= 0;
      // first order derivative approximation 
      fdm::Scalar h = mesh[1] - mesh[0];  
      Mat.coeffRef(0,0)=  val_coeff + deriv_coeff*(-1.0/h);
      Mat.coeffRef(0,1)=  deriv_coeff*(1.0/h);
    }; 
    void SetStencilR(fdm::Real t, const fdm::Vector& mesh, fdm::CSRMatrix& Mat) const 
    {
      Mat.bottomRows(1) *= 0; 
      // first order derivative approximation 
      fdm::Scalar h = mesh[mesh.size()-1] - mesh[mesh.size()-2];  
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-2)= deriv_coeff*(-1.0/h);
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=  val_coeff + deriv_coeff*(1.0/h);
    };

    void SetImpSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[0] = boundary_target;};
    void SetImpSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    {Sol[Sol.size()-1] = boundary_target;};
    
    // change the first/last (left/right boundary) entry of a vector  
    void SetSolL(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const 
    { 
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights(mesh[0], mesh.cbegin(), mesh.cbegin()+3, 1); 

      // solve the equation target = a*(S[0]) + b * ( W[0]*S[0] + W[1]*S[1] + W[2]*S[2] ) 
      // for the target value S[0] 
      double target = boundary_target; 
      target -= deriv_coeff * weights[1]*Sol[1]; 
      target -= deriv_coeff * weights[2]*Sol[2];
      target /= val_coeff + deriv_coeff * weights[0]; 

      // assign to Sol reference
      Sol[0] = target;  
      // void return type
    };
    void SetSolR(fdm::Real t, const fdm::Vector& mesh, fdm::StridedRef Sol) const  
    {
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      fdm::utils::FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights(mesh[mesh.size()-1], mesh.cend()-3, mesh.cend(), 1); 

      // solve the equation target = a*(S[N-1]) + b * ( W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] ) 
      // for the target value S[N-1] 
      fdm::Scalar target = boundary_target; 
      target -= deriv_coeff * weights[0]*Sol[Sol.size()-3]; 
      target -= deriv_coeff * weights[1]*Sol[Sol.size()-2]; 
      target /= val_coeff + deriv_coeff*weights[2]; 

      // assign to Sol reference
      Sol[Sol.size()-1] = target;  
      // void return type
    };
};

  } // end namespace osteps
} // end namespace fdm 

#endif // RobinBC.hpp