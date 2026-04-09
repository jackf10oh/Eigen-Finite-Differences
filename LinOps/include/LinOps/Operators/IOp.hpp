// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include<variant> 

#include "../LinOpMixin.hpp"
#include "../LinOpBase.hpp"

namespace linops{

class IOp : public LinOpMixIn<IOp>, public LinOpBase1D<IOp>, public LinOpBaseXD<IOp>
{
  private:
    // Member Data ---------------------------- 
    std::size_t m_size;
  public: 
    // Constructors + Destructor --------------------------
    IOp(std::size_t s_init=0) : m_size(s_init){} 
    IOp(const Mesh1D_SPtr_t& m){ setMesh1D(m); } 
    IOp(const MeshXD_SPtr_t& m){ setMeshXD(m);} 

    // destructor 
    ~IOp()=default;

    // Member Funcs  ======================================================

    // matrix getters 
    auto asMatrix() const { return Eigen::MatrixXd::Identity(m_size,m_size).sparseView(); }

    // Identity just returns inputs as outputs
    linops::Vector1D apply(const linops::Vector1D& d_arr) const { return d_arr; } 
    linops::VectorXD apply(const linops::VectorXD& d_arr) const { return d_arr; } 

    // fit operator to a domain mesh 
    void setMesh1D_impl(const Mesh1D_SPtr_t& m) 
    {
      m_size = m->size(); 
    };

    // fit operator to a domain mesh 
    void setMeshXD_impl(const MeshXD_SPtr_t& m) 
    {
      m_size = m->sizes_product(); 
    };

}; // End IOp 

} // end namespace linops 

#endif