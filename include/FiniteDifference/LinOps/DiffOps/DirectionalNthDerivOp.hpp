// DirectionalNthDerivOp.hpp
//
// X Dimensional version of NthDerivOp.hpp
//
// JAF 1/10/2026 

#ifndef DIRECTIONALNTHDERIVOP_H
#define DIRECTIONALNTHDERIVOP_H 

#include "../../Utilities/BlockDiagExpr.hpp"
#include "../../Utilities/HighDimExpr.hpp" 
#include "../LinOpBase.hpp" 
#include "NthDerivOp.hpp" 

namespace fdm{
  namespace linops{

template<std::size_t dir, std::size_t orderN>
class DirectionalNthDerivOp : public LinOpMixIn<DirectionalNthDerivOp<dir,orderN>>, public LinOpBaseXD<DirectionalNthDerivOp<dir,orderN>>
{
  private:
    // Member Data ----------------------------------------------
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 
    NthDerivOp<orderN> m_onedim_stencil;

  public:
    static constexpr std::size_t direction = dir; 
    static constexpr std::size_t order = orderN; 

    // Constructors + Destructor =====================================================
    // default
    DirectionalNthDerivOp()=default; 

    // meshxd 
    DirectionalNthDerivOp(const fdm::SharedConstMeshXD& m){this->setMesh(m);}

    // copy 
    DirectionalNthDerivOp(const linops::DirectionalNthDerivOp<dir,orderN>& other)=default; 

    // destructor 
    ~DirectionalNthDerivOp()=default; 

    // Member Functions =====================================================

    // Getters to Matrix 
    auto asMatrix() const { return make_HighDim(make_BlockDiag(  m_onedim_stencil.asMatrix(),m_prod_before),m_prod_after); }; 

  protected: 
    // Unreachable ------------------------------------------------------------
    // set DIrectional Derivative to operate on a domain mesh 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m)
    {
      m_onedim_stencil.set_mesh(m->getMesh1DSafe(dir)); // checks the safety that # of numDims >= Direction    
      m_prod_before = m->sizesMiddleProduct(dir+1,m->numDims()); // checks that m_dir < m->numDims()
      m_prod_after = m->sizesMiddleProduct(0,dir); 
    }
    
}; // end class DirectionalNthDerivOp

  } // end namespace linops 
} // end namespace fdm  

#endif // DirectionalNthDerivOp.hpp 

