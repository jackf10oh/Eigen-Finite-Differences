// DirectionalNthDerivOp.hpp
//
// X Dimensional version of NthDerivOp.hpp
//
// JAF 1/10/2026 

#ifndef DIRECTIONALNTHDERIVOP_H
#define DIRECTIONALNTHDERIVOP_H 

#include<Utilities/BlockDiagExpr.hpp> 
#include<Utilities/HighDimExpr.hpp> 

#include "../../LinearOpBase.hpp" 
#include "NthDerivOp.hpp" 

namespace linops{

class DirectionalNthDerivOp : public LinOpMixIn<DirectionalNthDerivOp>, public LinOpBaseXD<DirectionalNthDerivOp>
{
  private:
    // Member Data ----------------------------------------------
    MeshXD_WPtr_t m_mesh_ptr; 
    std::size_t m_dir; 
    std::size_t m_order; 
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 
    NthDerivOp m_onedim_stencil;

  public:
    // Constructors + Destructor =====================================================
    // meshxd + order + direction 
    DirectionalNthDerivOp(const MeshXD_SPtr_t& m, std::size_t order=1, std::size_t dir=0)
      : m_order(order), m_dir(dir), m_prod_before(std::size_t{1}), m_prod_after(std::size_t{1}), m_onedim_stencil(order) 
    {set_mesh(m);};

    // order + direction 
    DirectionalNthDerivOp(std::size_t order=1, std::size_t dir=0)
      : m_order(order), m_dir(dir), m_prod_before(1), m_prod_after(1), m_onedim_stencil(order) 
    {};

    // destructor 
    ~DirectionalNthDerivOp()=default; 

    // Member Functions =====================================================
    // getters to Direction + Order 
    std::size_t Direction() const {return m_dir; }
    std::size_t Order() const {return m_order; }

    // Getters to Matrix 
    auto asMatrix(){ return make_HighDim(make_BlockDiag( m_onedim_stencil.asMatrix(),m_prod_before),m_prod_after); }; 
    auto asMatrix() const { return make_HighDim(make_BlockDiag(  m_onedim_stencil.asMatrix(),m_prod_before),m_prod_after); }; 

    // getters to mesh 
    MeshXD_WPtr_t get_weak_meshxd() const { return m_mesh_ptr; }
    MeshXD_SPtr_t get_meshxd() const { return m_mesh_ptr.lock(); }

    // set DIrectional Derivative to operate on a domain mesh 
    void set_mesh(const MeshXD_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      if(!m) throw std::runtime_error("DirectionalNthDerivOp::set_mesh() error: shared_ptr<> expired!"); 
      // store the mesh   
      m_mesh_ptr = m; 

      // perform work on m ... 
      m_prod_before = m->sizes_middle_product(m_dir+1,m->dims()); // checks that m_dir < m->dims()
      m_prod_after = m->sizes_middle_product(0,m_dir); 
      m_onedim_stencil.set_mesh(m->GetMesh(m_dir)); // no need for range safe ->GetMeshAt()  
    }
    
}; // end class DirectionalNthDerivOp

} // end namespace linops 

#endif // DirectionalNthDerivOp.hpp 

