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
    AutonomousCoeff(const Mesh1D_SPtr_t m, Callable f)
      : CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMesh1D(m);}

    // from Mesh1D + callable 
    AutonomousCoeff(const MeshXD_SPtr_t m, Callable f)
      : CoeffOpMixIn<AutonomousCoeff<Callable>, Callable>(f) 
    {this->setMeshXD(m);}

    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=default; 
    
    // destructor
    ~AutonomousCoeff()=default;

    // Member Funcs =============================================================
    void setMesh1D_impl(const Mesh1D_SPtr_t& m)
    {
      CoeffOpMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }

    void setMeshXD_impl(const MeshXD_SPtr_t& m)
    {
      CoeffOpMixIn<AutonomousCoeff,Callable>::fillDiagonal(m); 
    }
}; 

} // end namespace linops 

#endif // AutonomousCoeff.hpp 