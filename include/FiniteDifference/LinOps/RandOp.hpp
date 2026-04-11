// RandOp.hpp
//
// header file for a linop wrapper of a random matrix
//
// JAF 12/6/2025

#ifndef RANDLINOP_H
#define RANDLINOP_H

#include<Eigen/Core>
#include "LinOpBase.hpp"
#include "LinOpMixin.hpp" 

namespace linops{

class RandOp : public LinOpMixIn<RandOp>, public LinOpBase1D<RandOp>, public LinOpBaseXD<RandOp>
{
  private:
    // Member Data -------------------------------  
    fdm::Matrix m_mat; 
  public: 

    // Constructors + Destructor ===========================
    RandOp()=default;
    RandOp(const fdm::SharedConstMesh1D& m){ setMesh(m); }   
    RandOp(const fdm::SharedConstMeshXD& m){ setMesh(m); }
    RandOp(const RandOp& other)=default; 
    
    // destructor
    ~RandOp()=default;

    // Member Funcs ==========================================
    using LinOpBase1D<RandOp>::setMesh; 
    using LinOpBase1D<RandOp>::apply; 
    using LinOpBaseXD<RandOp>::setMesh; 
    using LinOpBaseXD<RandOp>::apply; 

    // matrix getters 
    const auto& asMatrix() const { return m_mat; };
    
  protected: 
    // Unreachable ------------------------------------------------------------
    // fit operator to a domain mesh 
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m)
    {
      auto s = m->size();       
      m_mat = Eigen::MatrixXd::Random(s, s).sparseView(); 
    };

    // fit operator to a domain mesh 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m)
    {
      auto s = m->sizesProduct(); 
      m_mat = Eigen::MatrixXd::Random(s, s).sparseView(); 
    };

}; // end  RandOp

} // end namespace linops 

#endif // RandOp.hpp