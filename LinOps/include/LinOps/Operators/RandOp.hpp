// RandOp.hpp
//
// header file for a linop wrapper of a random matrix
//
// JAF 12/6/2025

#ifndef RANDLINOP_H
#define RANDLINOP_H

#include<Eigen/Core>

#include "../LinOpBase.hpp"
#include "../LinOpMixin.hpp" 

namespace linops{

class RandOp : public LinOpMixIn<RandOp>, public LinOpBase1D<RandOp>, public LinOpBaseXD<RandOp>
{
  private:
    // Member Data -------------------------------  
    linops::Matrix m_mat; 
  public: 

    // Constructors + Destructor ===========================
    RandOp(const Mesh1D_SPtr_t& m){ setMesh1D(m); }   
    RandOp(const MeshXD_SPtr_t& m){ setMeshXD(m); }
    
    // destructor
    ~RandOp()=default;

    // Member Funcs ==========================================

    // matrix getters 
    auto asMatrix() const { return m_mat; };

    // fit operator to a domain mesh 
    void setMesh1D_impl(const Mesh1D_SPtr_t& m)
    {
      m_mat = Eigen::MatrixXd::Random(m->size(), m->size()).sparseView(); 
    };

    // fit operator to a domain mesh 
    void setMeshXD_impl(const MeshXD_SPtr_t& m)
    {
      m_mat = Eigen::MatrixXd::Random(m->sizes_product(), m->sizes_product()).sparseView(); 
    };

}; // end  RandOp

} // end namespace linops 

#endif // RandOp.hpp