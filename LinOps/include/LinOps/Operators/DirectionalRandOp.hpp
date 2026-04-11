// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H

#include<Utilities/BlockDiagExpr.hpp> 
#include<Utilities/HighDimExpr.hpp> 

#include "../LinOpMixIn.hpp"
#include "../LinOpBase.hpp" 

namespace linops{

class DirectionalRandOp: public LinOpMixIn<DirectionalRandOp>, public LinOpBaseXD<DirectionalRandOp> 
{
  private:
    private:
      // Member Data ---------------------------
      std::size_t m_direction; // which Mesh1D the operator acts on. 
      std::size_t m_prod_before; // which Mesh1D the operator acts on. 
      std::size_t m_prod_after; // which Mesh1D the operator acts on. 
      linops::Matrix m_mat; 
  public:
    // Constructors + Destructo =====================================

    // direction only 
    DirectionalRandOp(std::size_t dir_init=0) 
      : m_direction(dir_init)
    {};

    // mesh + direction 
    DirectionalRandOp(std::size_t dir_init, const SharedConstMeshXD& m) 
      : m_direction(dir_init) 
    {setMesh(m);}; 

    // destructor
    ~DirectionalRandOp()=default; 

    // Member Funcs ============================================== 
    
    // matrix getters
    auto asMatrix() const { return make_HighDim(make_BlockDiag(m_mat,m_prod_before),m_prod_after); }; 

  protected: 
    // Unreachable ------------------------------------------------------------
    // set operator to domain mesh 
    void setMeshXD_impl(const SharedConstMeshXD& m){
      m_prod_before = m->sizesMiddleProduct(m_direction+1,m->numDims()); 
      m_prod_after = m->sizesMiddleProduct(0,m_direction); 
      m_mat =  Eigen::MatrixXd::Random(m->sizeOfDim(m_direction),m->sizeOfDim(m_direction)).sparseView(); 
    }  
}; 

} // end namespace linops 

#endif 