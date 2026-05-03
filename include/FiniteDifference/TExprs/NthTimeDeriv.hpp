// NthTimeDeriv.hpp
//
//
//
// JAF 1/15/2026 

#ifndef NTHTIMEDERIV_H
#define NTHTIMEDERIV_H 

#include "TimeDerivBase.hpp"

namespace fdm{
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
} // end namespace fdm 

#endif // NthTimeDeriv.hpp