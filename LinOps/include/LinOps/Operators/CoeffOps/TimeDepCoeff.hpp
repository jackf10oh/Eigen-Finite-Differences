// TimeDepCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include "../../CoeffOpMixIn.hpp" 
#include "../../LinOpBase.hpp" 
#include "../../LinOpTraits.hpp" // callable traits -> bind_first

namespace linops{

template<typename Callable>
class TimeDepCoeff : public CoeffOpMixIn<TimeDepCoeff<Callable>, typename traits::callable_traits<Callable>::BindFirst_t>, public LinOpBase1D<TimeDepCoeff<Callable>>, public LinOpBaseXD<TimeDepCoeff<Callable>>
{
  private:
    // convenience type
    using Binded = typename traits::callable_traits<Callable>::BindFirst_t; 
    // Member Data -------------------------- 
    linops::MeshXD_SPtr_t m_owned_mesh; 
  public:
    // Constructors + Destructor ==========================================================
    // no default constructor
    TimeDepCoeff()=delete; 

    // from callable 
    TimeDepCoeff(Callable f)
      : CoeffOpMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {}

    // from Mesh1D + callable 
    TimeDepCoeff(const Mesh1D_SPtr_t m, Callable f)
      : CoeffOpMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {this->setMesh1D(m);}

    // from Mesh1D + callable 
    TimeDepCoeff(const MeshXD_SPtr_t m, Callable f)
      : CoeffOpMixIn<TimeDepCoeff<Callable>, Binded>(Binded(f, 0.0)) 
    {this->setMeshXD(m);}

    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other)=default; 
    
    // destructor
    ~TimeDepCoeff()=default;

    // Member Funcs =============================================================
    void setMesh1D_impl(const Mesh1D_SPtr_t& m)
    {
      m_owned_mesh = linops::make_MeshXD(m); 
    }

    void setMeshXD_impl(const MeshXD_SPtr_t& m)
    {
      m_owned_mesh = linops::make_MeshXD(m); 
    }

    void setTime_impl(double t)
    {
      this->m_callable.captured = t; 
      this->fillDiagonal(m_owned_mesh); 
    }
}; 

} // end namespace linops 

#endif // TimeDepCoeff.hpp 