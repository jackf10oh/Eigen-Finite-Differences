// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H 

#include "LinOpMixIn.hpp"
#include "LinOpBase.hpp" 
#include "../Utilities/BlockDiagExpr.hpp"
#include "../Utilities/HighDimExpr.hpp"

namespace linops{

class DirectionalRandOp: public LinOpMixIn<DirectionalRandOp>, public LinOpBaseXD<DirectionalRandOp> 
{
  private:
    private:
      // Member Data ---------------------------
      std::size_t m_direction; // which Mesh1D the operator acts on. 
      std::size_t m_prod_before; // which Mesh1D the operator acts on. 
      std::size_t m_prod_after; // which Mesh1D the operator acts on. 
      fdm::Matrix m_mat; 
  public:
    // Constructors + Destructo =====================================

    // direction only 
    DirectionalRandOp(std::size_t dir_init=0) 
      : m_direction(dir_init)
    {};

    // mesh + direction 
    DirectionalRandOp(std::size_t dir_init, const fdm::SharedConstMeshXD& m) 
      : m_direction(dir_init) 
    {setMesh(m);}; 

    // destructor
    ~DirectionalRandOp()=default; 

    // Member Funcs ============================================== 
    
    // matrix getters
    auto asMatrix() const { return fdm::utils::make_HighDim(fdm::utils::make_BlockDiag(m_mat,m_prod_before),m_prod_after); }; 

  protected: 
    // Unreachable ------------------------------------------------------------
    // set operator to domain mesh 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m){
      m_prod_before = m->sizesMiddleProduct(m_direction+1,m->numDims()); 
      m_prod_after = m->sizesMiddleProduct(0,m_direction); 
      m_mat =  Eigen::MatrixXd::Random(m->sizeOfDim(m_direction),m->sizeOfDim(m_direction)).sparseView(); 
    }  
}; 

} // end namespace linops 

#endif 