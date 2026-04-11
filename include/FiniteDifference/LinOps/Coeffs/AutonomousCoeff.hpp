// AutonomousCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include "../CoeffMixIn.hpp" 
#include "../LinOpBase.hpp" 

namespace linops{

template<typename Callable>
class AutonomousCoeff : public CoeffMixIn<AutonomousCoeff<Callable>, Callable>, public LinOpBase1D<AutonomousCoeff<Callable>>, public LinOpBaseXD<AutonomousCoeff<Callable>>
{
  public:
    // Constructors + Destructor ==========================================================
    // no default constructor
    AutonomousCoeff()=delete; 

    // from callable 
    AutonomousCoeff(Callable f)
      : CoeffMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {}

    // from Mesh1D + callable 
    AutonomousCoeff(Callable f, const fdm::SharedConstMesh1D& m)
      : CoeffMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMesh(m);}

    // from Mesh1D + callable 
    AutonomousCoeff(Callable f, const fdm::SharedConstMeshXD& m)
      : CoeffMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMesh(m);}

    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=default; 
    
    // destructor
    ~AutonomousCoeff()=default;

    // Member Funcs =============================================================
    using LinOpBase1D<AutonomousCoeff>::setMesh; 
    using LinOpBase1D<AutonomousCoeff>::apply; 
    using LinOpBaseXD<AutonomousCoeff>::setMesh; 
    using LinOpBaseXD<AutonomousCoeff>::apply; 

  protected: 
    // Unreachable ------------------------------------------------------------
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m)
    {
      CoeffMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }

    void setMeshXD_impl(const fdm::SharedConstMeshXD& m)
    {
      CoeffMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }
}; 

} // end namespace linops 

#endif // AutonomousCoeff.hpp 