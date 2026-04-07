// NthTimeDeriv.hpp
//
//
//
// JAF 1/15/2026 

#ifndef NTHTIMEDERIV_H
#define NTHTIMEDERIV_H 

#include "TimeDerivBase.hpp"

namespace TExprs{

// ===============================================================
template<std::size_t N>
class NthTimeDeriv : public TExprs::TimeDerivBase<NthTimeDeriv<N>, N> 
{

  public:
    // Memeber Data ------------------------------------- 
    static constexpr std::size_t order = N; 

    // Constructors + Destructor ====================================
    NthTimeDeriv()=default;
    NthTimeDeriv(const NthTimeDeriv& other)=default;
    
    // destructor 
    ~NthTimeDeriv()=default; 

    // Member Funcs =================================================== 
    template<std::size_t ithCol, std::size_t nCols, typename Cont>
    decltype(auto) coeffAt(const Cont& v) const 
    {
      return v[nCols * order + ithCol]; 
    }

    // using LhsBase<NthTimeDeriv>::toTuple; 
    std::string toString() const {return "hi from NthTimeDeriv"; }; 

}; 

} // end namespace TExprs 

#endif // NthTimeDeriv.hpp