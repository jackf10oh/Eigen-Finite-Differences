// Robin.hpp 
//
// robin(third type) boundary conditions 
// alpha*U + beta*Ux = target_val
//
// JAF 12/8/2025

#ifndef FORNFDM_OSTEPS_ROBINBC_H
#define FORNFDM_OSTEPS_ROBINBC_H

#include "../../Utilities/FornbergStackCalc.hpp"
#include "BCPair.hpp"

namespace fornfdm{
  namespace osteps{

class RobinBC
{
  public:  
    // member data 
    fornfdm::Scalar flux;
    fornfdm::Scalar stiffness;
    fornfdm::Scalar damping;

  public:
    // Constructors + Destructor ---------------------------------------------
    // a*U + b*Ux = target
    RobinBC(fornfdm::Scalar a=1.0, fornfdm::Scalar b=0.0, fornfdm::Scalar target=0.0) 
      : stiffness(a), damping(b), flux(target)
    {}; 
    // copy 
    RobinBC(const RobinBC& other)=default; 
    // destructors
    virtual ~RobinBC()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fornfdm stencil matrix
    BoundaryRow getTopRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    {
      // first order derivative approximation 
      fornfdm::Scalar h = axis[1] - axis[0];  
      return BoundaryRow{{stiffness + damping*(-1.0/h), damping*(1.0/h), 0},2};
    }; 
    BoundaryRow getBottomRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    {
      // first order derivative approximation 
      fornfdm::Scalar h = axis[1] - axis[0];  
      return BoundaryRow{{damping*(-1.0/h), stiffness + damping*(1.0/h), 0}, 2};
    }; 

    void setImpSolLeft(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef sol) const 
    {sol[0] = flux;};
    void setImpSolRight(fornfdm::Real t, const fornfdm::Vector& mesh, fornfdm::StrideRef sol) const 
    {sol[sol.size()-1] = flux;};
    
    // change the first/last (left/right boundary) entry of a vector  
    void setExpSolLeft(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const 
    { 
      // up to 3 nodes, up to 1st order deriv
      fornfdm::utils::FornbergStackCalc<3,1> calc;
      calc.calculate(axis[0], axis.cbegin(), std::next(axis.cbegin(),3));

      // get forward finite difference weights for sol[0], sol[1], sol[2] 
      // solve the equation target = a*(S[0]) + b * ( W[0]*S[0] + W[1]*S[1] + W[2]*S[2] ) 
      // for the target value S[0] 
      double target = flux; 
      target -= damping * calc.getArray()[4]*sol[1]; 
      target -= damping * calc.getArray()[5]*sol[2];
      target /= stiffness + damping * calc.getArray()[3]; 

      // assign to sol reference
      sol[0] = target;  
    };
    void setExpSolRight(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const  
    {
      // up to 3 nodes, up to 1st order deriv
      fornfdm::utils::FornbergStackCalc<3,1> calc;
      calc.calculate(axis[axis.size()-1], std::prev(axis.cend(),3), axis.cend());

      // solve the equation target = a*(S[N-1]) + b * ( W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] ) 
      // for the target value S[N-1] 
      fornfdm::Scalar target = flux; 
      target -= damping * calc.getArray()[3] * sol[sol.size()-3]; 
      target -= damping * calc.getArray()[4] * sol[sol.size()-2]; 
      target /= stiffness + damping * calc.getArray()[5]; 

      // assign to sol reference
      sol[sol.size()-1] = target;  
      // void return type
    };
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // RobinBC.hpp