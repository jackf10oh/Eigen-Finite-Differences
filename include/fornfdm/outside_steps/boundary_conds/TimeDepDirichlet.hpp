// TimeDepDirichlet.hpp 
//
// dirichlet boundary value conditions of the form 
// U(0) = c(t) for some callable c 
//
// JAF 6/6/2026

#ifndef FORNFDM_OSTEPS_TIMEDEPDIRICHLET_H
#define FORNFDM_OSTEPS_TIMEDEPDIRICHLET_H

#include "BCPair.hpp" // struct BoundaryRow

namespace fornfdm{
  namespace osteps{

template<typename Callable>
class TimeDepDirichlet
{
  public:  
    // member data 
    std::remove_reference_t<Callable> callable;

  public:
    // Constructors ---------------------------------------------
    TimeDepDirichlet(Callable c) : callable(c){} 
    TimeDepDirichlet(const TimeDepDirichlet& other)=default; 
    ~TimeDepDirichlet()=default; 

    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fornfdm stencil matrix
    BoundaryRow getTopRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    { return BoundaryRow{ {1,0,0}, 1}; }

    BoundaryRow getBottomRow(fornfdm::Real t, const fornfdm::Vector& axis) const 
    { return BoundaryRow{ {1,0,0}, 1}; }

    // change the first/last (left/right boundary) entry of a vector to implicit solution   
    void setImpSolLeft(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const 
    {sol[0] = callable(t); }

    void setImpSolRight(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const 
    {sol[sol.size()-1] = callable(t); }
    
    // change the first/last (left/right boundary) entry of a vector to explicit solution
    void setExpSolLeft(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const 
    { sol[0] = callable(t); }

    void setExpSolRight(fornfdm::Real t, const fornfdm::Vector& axis, fornfdm::StrideRef sol) const 
    {sol[sol.size()-1] = callable(t); }
};

  } // end namespace osteps
} // end namespace fornfdm 

#endif // TimeDepDirichlet.hpp