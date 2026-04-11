// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include "../Utilities/RowMajorIdentityExpr.hpp"
#include "LinOpMixin.hpp"
#include "LinOpBase.hpp"

namespace linops{

class IOp : public LinOpMixIn<IOp>, public LinOpBase1D<IOp>, public LinOpBaseXD<IOp>
{
  private:
    // Member Data ---------------------------- 
    std::size_t m_size;     
  public: 
    // Constructors + Destructor --------------------------
    IOp(std::size_t s_init=0) : m_size(s_init){} 
    IOp(const fdm::SharedConstMesh1D& m){ setMesh(m); } 
    IOp(const fdm::SharedConstMeshXD& m){ setMesh(m);} 

    // destructor 
    ~IOp()=default;

    // Member Funcs  ======================================================
    using LinOpBase1D<IOp>::setMesh; 
    using LinOpBase1D<IOp>::apply; 
    using LinOpBaseXD<IOp>::setMesh; 
    using LinOpBaseXD<IOp>::apply; 

    // matrix getters 
    auto asMatrix() const { return fdm::utils::make_RowMajorIdentity(m_size,m_size); }; 

    // Identity just returns inputs as outputs
    fdm::Vector1D apply(const fdm::Vector1D& d_arr) const { return d_arr; } 
    fdm::VectorXD apply(const fdm::VectorXD& d_arr) const { return d_arr; } 

  protected: 
    // Unreachable ------------------------------------------------------------
    // fit operator to a domain mesh 
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m) 
    {
      m_size = m->size(); 
    };

    // fit operator to a domain mesh 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m) 
    {
      m_size = m->sizesProduct();
    };

}; // End IOp 

} // end namespace linops 

#endif