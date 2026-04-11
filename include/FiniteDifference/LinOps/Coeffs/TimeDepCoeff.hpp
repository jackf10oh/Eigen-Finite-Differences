// TimeDepCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include "../CoeffMixIn.hpp" 
#include "../LinOpBase.hpp" 
#include "../LinOpTraits.hpp" // callable traits -> bind_first

namespace fdm{
  namespace linops{

template<typename Callable>
class TimeDepCoeff : public CoeffMixIn<TimeDepCoeff<Callable>, typename traits::callable_traits<Callable>::BindFirst_t>, public LinOpBase1D<TimeDepCoeff<Callable>>, public LinOpBaseXD<TimeDepCoeff<Callable>>
{
  private:
    // convenience type
    using Binded = typename traits::callable_traits<Callable>::BindFirst_t; 
    // Member Data -------------------------- 
    fdm::SharedConstMeshXD m_owned_mesh; 
  public:
    // Constructors + Destructor ==========================================================
    // no default constructor
    TimeDepCoeff()=delete; 

    // from callable 
    TimeDepCoeff(Callable f)
      : CoeffMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {}

    // from Mesh1D + callable 
    TimeDepCoeff(Callable f, const fdm::SharedConstMesh1D& m)
      : CoeffMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {this->setMesh(m);}

    // from Mesh1D + callable 
    TimeDepCoeff(Callable f, const fdm::SharedConstMeshXD& m)
      : CoeffMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {this->setMesh(m);}

    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other)=default; 
    
    // destructor
    ~TimeDepCoeff()=default;

    // Member Funcs =============================================================
    using LinOpBase1D<TimeDepCoeff>::setMesh; 
    using LinOpBase1D<TimeDepCoeff>::apply; 
    using LinOpBaseXD<TimeDepCoeff>::setMesh; 
    using LinOpBaseXD<TimeDepCoeff>::apply; 
    
  protected: 
    // Unreachable ------------------------------------------------------------
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m)
    {
      m_owned_mesh = fdm::make_MeshXD(m); 
    }

    void setMeshXD_impl(const fdm::SharedConstMeshXD& m)
    {
      m_owned_mesh = fdm::make_MeshXD(m); 
    }

    void setTime_impl(double t)
    {
      this->m_callable.captured = t; 
      this->fillDiagonal(m_owned_mesh); 
    }
}; 

  } // end namespace linops 
} // end namespace fdm  

#endif // TimeDepCoeff.hpp 