// AutonomousCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include "../../CoeffOpMixIn.hpp" 
#include "../../LinOpBase.hpp" 

namespace linops{

template<typename Callable>
class AutonomousCoeff : public CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>, public LinOpBase1D<AutonomousCoeff<Callable>>, public LinOpBaseXD<AutonomousCoeff<Callable>>
{
  public:
    // Constructors + Destructor ==========================================================
    // no default constructor
    AutonomousCoeff()=delete; 

    // from callable 
    AutonomousCoeff(Callable f)
      : CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {}

    // from Mesh1D + callable 
    AutonomousCoeff(const SharedConstMesh1D m, Callable f)
      : CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMesh1D(m);}

    // from Mesh1D + callable 
    AutonomousCoeff(const SharedConstMeshXD m, Callable f)
      : CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMesh(m);}

    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=default; 
    
    // destructor
    ~AutonomousCoeff()=default;

    // Member Funcs =============================================================
    using LinOpBase1D<AutonomousCoeff>::setMesh1D; 
    using LinOpBase1D<AutonomousCoeff>::apply; 
    using LinOpBaseXD<AutonomousCoeff>::setMesh; 
    using LinOpBaseXD<AutonomousCoeff>::apply; 

    void setMesh1D_impl(const SharedConstMesh1D& m)
    {
      CoeffOpMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }

    void setMeshXD_impl(const SharedConstMeshXD& m)
    {
      CoeffOpMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }
}; 

} // end namespace linops 

#endif // AutonomousCoeff.hpp 