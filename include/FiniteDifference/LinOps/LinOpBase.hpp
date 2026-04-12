// LineOpBase.hpp
//
// CTRP base for 1D differential operator L
// where L works on function discretizations across x0 , ..., xN 
//
// JAF 12/7/2025

#ifndef LINEAROPBASE_H
#define LINEAROPBASE_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>

#include "LinOpTraits.hpp"

#include "../Mesh1D.hpp"
#include "../Vector1D.hpp"
#include "../MeshXD.hpp"
#include "../VectorXD.hpp"

namespace fdm{
  namespace linops{

// forward declaration 
template<typename T> class LinOpBaseXD; 

template<typename Derived>
class LinOpBase1D 
{
  friend LinOpBaseXD<Derived>; 
  public:
    // Type Defs -----------------------------------------------------
    typedef struct{} is_1dim_linop_tag; // to tell if a class derived from LinOpBase1D<>
  
  private:
    // Member Data ---------------------------------------------------- 
    fdm::WeakConstMesh1D m_mesh_ptr; 

  public:
    // Member Funcs  ======================================================

    // multiply the underlying expression with linops::Vector1D's underlying Eigen::VectorXd
    fdm::Vector1D apply(const fdm::Vector1D& d) const 
    {
      Eigen::VectorXd v = static_cast<const Derived*>(this)->asMatrix() * d.values();  // calculate A*b
      fdm::Vector1D result(std::move(v), this->m_mesh_ptr); // move A*b into result's values
      return result;
    };

    // fit operator to a mesh of rectangular domain 
    void setMesh(const fdm::SharedConstMesh1D& m) 
    {
      // if Derived is an expression 
      if constexpr(traits::is_expr_crtp<Derived>::value){
        auto& expr = static_cast<Derived&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename Derived::LStorage>::value) expr.getLhs().setMesh(m); 
        if constexpr(traits::is_1dim_linop_crtp<typename Derived::RStorage>::value) expr.getRhs().setMesh(m); 
      }
      // Non Expression case ... 
      else{
        /* // // ensure we aren't resetting the mesh again
        // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
        // // on nullptr throw an error  
        // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
        */ 
       // store the mesh
        m_mesh_ptr = m;   
        // if the derived is also a 1D linop clear the mesh 
        if constexpr(traits::is_xdim_linop_crtp<Derived>::value){
          static_cast<Derived*>(this)->clearMeshXD();
        }
        // // perform work on m 
        setMesh1D_impl(m);
      }
    };
    
    // return Mesh1D pointed to
    fdm::SharedConstMesh1D getMesh1D() const 
    {
      // if Derived is an expression 
      if constexpr(traits::is_expr_crtp<Derived>::value){
        auto& expr = static_cast<const Derived&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename Derived::LStorage_t>::value) return expr.getLhs().getMesh1D(); 
        else if constexpr(traits::is_1dim_linop_crtp<typename Derived::RStorage_t>::value) return expr.getRhs().getMesh1D(); 
        else static_assert(false, "cannot call getMesh1D() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_mesh_ptr.lock();
      }
    } 

  protected:
    // Unreachable except by other base. 
    void clearMesh1D(){ this->m_mesh_ptr = fdm::SharedConstMesh1D{}; }

    // Must implement! -------------------------
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m)
    {
      // static_cast<Derived*>(this)->setMesh1D_impl(m);
      Accessor::call_setMesh1D_impl(*this, m);  
    }

    struct Accessor : public Derived
    {
      static void call_setMesh1D_impl(LinOpBase1D& self, const fdm::SharedConstMesh1D& m){ static_cast<Accessor&>(self).setMesh1D_impl(m); }
    }; 

}; // end LinOpBase1D<>  

template<typename Derived>
class LinOpBaseXD 
{
  friend LinOpBase1D<Derived>; 
  public:
    // Type Defs -----------------------------------------------------
    typedef struct{} is_xdim_linop_tag; // to tell if a class derived from LinOpBase1D<>

  private:
    // Member Data ------------------------------------------------- 
    fdm::WeakConstMeshXD m_mesh_ptr; 

  public:
    // Member Funcs  ======================================================

    // multiply the underlying expression with linops::VectorXD's underlying Eigen::VectorXd
    fdm::VectorXD apply(const fdm::VectorXD& d) const 
    {
      Eigen::VectorXd v = static_cast<const Derived*>(this)->asMatrix() * d.values();  // calculate A*b
      fdm::VectorXD result(std::move(v), this->m_mesh_ptr); // move A*b into result's values
      return result;
    };

    // fit operator to a mesh of rectangular domain 
    void setMesh(const fdm::SharedConstMeshXD& m) 
    {
      // if Derived is an expression 
      if constexpr(traits::is_expr_crtp<Derived>::value){
        auto& expr = static_cast<Derived&>(*this);  
        if constexpr(traits::is_xdim_linop_crtp<typename Derived::LStorage_t>::value) expr.getLhs().setMesh(m); 
        if constexpr(traits::is_xdim_linop_crtp<typename Derived::RStorage_t>::value) expr.getRhs().setMesh(m); 
      }
      // Non Expression case ... 
      else{
        /* // // ensure we aren't resetting the mesh again
        // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
        // // on nullptr throw an error  
        // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
        */ 
        // store the mesh
        m_mesh_ptr = m;  
        // if the derived is also a 1D linop clear the mesh 
        if constexpr(traits::is_1dim_linop_crtp<Derived>::value){
          static_cast<Derived*>(this)->clearMesh1D();
        }  
        // // perform work on m 
        // static_cast<Derived*>(this)->setMeshXD_impl(m);
        Accessor::call_setMeshXD_impl(*this, m); 
      }
    };
    
    // return Mesh1D pointed to
    fdm::SharedConstMeshXD getMeshXD() const 
    {
      // if Derived is an expression 
      if constexpr(traits::is_expr_crtp<Derived>::value){
        auto& expr = static_cast<const Derived&>(*this);  
        if constexpr(traits::is_xdim_linop_crtp<typename Derived::LStorage_t>::value) return expr.getLhs().getMeshXD(); 
        else if constexpr(traits::is_xdim_linop_crtp<typename Derived::RStorage_t>::value) return expr.getRhs().getMeshXD(); 
        else static_assert(false, "cannot call getMeshXD() on expr with no LinOpBaseXD's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_mesh_ptr.lock();
      }
    } 
    
  protected:
    // Unreachable except other base. 
    void clearMeshXD(){ this->m_mesh_ptr = fdm::SharedConstMeshXD{}; }

    // must implement! 
    void setMeshXD_impl(const fdm::SharedConstMeshXD& m)
    {
      // static_cast<Derived*>(this)->setMeshXD_impl(m); 
      Accessor::call_setMeshXD_impl(*this, m); 
    }

    struct Accessor : public Derived
    {
      static void call_setMeshXD_impl(LinOpBaseXD& self, const fdm::SharedConstMeshXD& m){ static_cast<Accessor&>(self).setMeshXD_impl(m); }
    }; 

}; // end LinOpBaseXD<>  

  } // end namespace linops 
} // end namespace fdm  

#endif // LinearOpBase.hpp