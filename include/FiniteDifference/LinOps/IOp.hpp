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

namespace fdm{
  namespace linops{

class IOp : public LinOpMixIn<IOp>, public LinOpBase1D<IOp>, public LinOpBaseXD<IOp>
{
  private:
    // Member Data ---------------------------- 
    std::size_t m_rows;    
    std::size_t m_cols;  
  public: 
    // Constructors + Destructor --------------------------
    IOp(std::size_t s_init=0) : m_rows(s_init), m_cols(s_init){} 
    IOp(std::size_t m, std::size_t n) : m_rows(m), m_cols(n){} 
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
    auto asMatrix() const { return fdm::utils::make_RowMajorIdentity(m_rows,m_cols); }; 

    // Identity just returns inputs as outputs
    fdm::Vector1D apply(const fdm::Vector1D& d_arr) const { return d_arr; } 
    fdm::VectorXD apply(const fdm::VectorXD& d_arr) const { return d_arr; } 

    // able to resize. useful for implicit solvers 
    void resize(std::size_t s){ m_rows = m_cols = s; }
    void resize(std::size_t m, std::size_t n){ m_rows = m; m_cols = n; }

  protected: 
    // Unreachable ------------------------------------------------------------
    // fit operator to a domain mesh 
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m) 
    {
      m_rows = m_cols = m->size(); 
    };

    // fit operator to a domain mesh 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m) 
    {
      m_rows = m_cols = m->sizesProduct();
    };

}; // End IOp 

  } // end namespace linops 
} // end namespace fdm  

#endif