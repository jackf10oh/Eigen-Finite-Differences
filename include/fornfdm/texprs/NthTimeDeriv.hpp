// NthTimeDeriv.hpp
//
// Concrete class derived from TimeDerivBase 
// That represent operator like Ut, Utt, etc 
// for time derivative template expressions 
//
// JAF 1/15/2026 

#ifndef FORNFDM_TEXPRS_NTHTIMEDERIV_H
#define FORNFDM_TEXPRS_NTHTIMEDERIV_H

#include<cstdint>
#include "traits.hpp"
#include "TimeDerivBase.hpp"

namespace fornfdm{
namespace texprs{

// ===============================================================
template<std::size_t N>
class NthTimeDeriv : public texprs::TimeDerivBase<NthTimeDeriv<N>, N> 
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
    auto coeffAt(const Cont& v) const 
    {
      return v[nCols * order + ithCol]; 
    }

}; 

} // end namespace texprs 
} // end namespace fornfdm 

#endif // NthTimeDeriv.hpp