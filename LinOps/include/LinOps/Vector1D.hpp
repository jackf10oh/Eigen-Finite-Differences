// Vector1D.hpp
//
// Discretization of a function F(x0), ..., F(xN) on a mesh x0, ..., xN
//
// JAF 12/5/2025

#ifndef VECTOR1D_H
#define VECTOR1D_H

#include<memory>
#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include <type_traits>

#include "Mesh1D.hpp"
#include "LinOpTraits.hpp"

namespace linops{

class Vector1D
{
  private:
    // Member Data --------------------------------------------------------------------
    Mesh1D_WPtr_t m_mesh_ptr; // weak pointer to mesh the function "maps" into the reals
    Eigen::VectorXd m_vals; // eigen array of functions vals at mesh points
  public:
    // Constructors / Destructors ==============================================================
    // Args -------------------------------
    // from size 
    Vector1D(std::size_t size_init=0): m_mesh_ptr(), m_vals(size_init){};
    // from mesh
    // from const shared_ptr<Mesh1D>
    Vector1D(const Mesh1D_SPtr_t& mesh_init) : m_mesh_ptr(mesh_init), m_vals(mesh_init->size()){}
    
    // Copy -----------------------------
    Vector1D(const Vector1D& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {};
    
    // copy from Eigen::VectorXd
    Vector1D(const Eigen::VectorXd& other, Mesh1D_WPtr_t mesh_init = Mesh1D_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(other)
    {}; 
    
    // Move ---------------------------
    Vector1D(Vector1D&& other)
      : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {}; 
    
    // move from Eigen::VectorXD
    Vector1D(Eigen::VectorXd&& other, Mesh1D_WPtr_t mesh_init = Mesh1D_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(std::move(other))
    {}; 
    // destructors ----------------------------------
    ~Vector1D()=default; 

    // Member Functions ==============================================================
    // return the size of the underlying vector ------------------------
    std::size_t size() const { return m_vals.size(); }

    // get underlying values ------------------
    Eigen::VectorXd& values(){ return m_vals; }
    const Eigen::VectorXd& values() const { return m_vals; }
    
    // indexing ------------------------------------
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ assert(i<size()); return m_vals[i]; }
    const double& at(std::size_t i) const { assert(i<size()); return m_vals[i]; }

    // get underlying mesh ------------------------------
    // const auto mesh(){return m_mesh_ptr.lock(); }; 
    Mesh1D_SPtr_t get_mesh1d() const{return m_mesh_ptr.lock(); }

    // set the mesh the discretization is on -----------
    auto& set_mesh(Mesh1D_WPtr_t m){ 
      m_mesh_ptr=m;
      return *this; 
    }
    // set vector to same size as mesh -----------------
    Vector1D& resize(const Mesh1D_SPtr_t& m){
      m_mesh_ptr=m;
      m_vals.conservativeResize(m->size()); 
      return *this; 
    }
    
    // forward iterators -------------------------------------------------------
    auto begin(){ return m_vals.begin(); }
    auto end(){ return m_vals.end(); }

    // forward citerators ----------------------------------
    auto cbegin() const { return m_vals.cbegin(); }
    auto cend() const { return m_vals.cend(); }

    // no reverse citerators in eigen 

    // Operators ================================================================
    // Assignment from Vector1D --------------------------------
    Vector1D& operator=(const Vector1D& other) = default;
    Vector1D& operator=(Vector1D&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      m_vals=std::move(other.m_vals); 
      return *this;
    }; 
}; // end Vector1D 

// set vector to a constant -------------------------------------------------------------
linops::Vector1D make_Discretization(const Mesh1D_SPtr_t& m, double val){
  linops::Vector1D result(m);
  result.values().setConstant(val); 
  return result; 
}

// resize + set vector F(x0),...,F(xN) for x0,...,xN in stored mesh 
template<
typename F,
typename = std::enable_if_t<
  !std::is_arithmetic_v<std::remove_reference_t<std::remove_cv_t<F>>>
  >
>
linops::Vector1D make_Discretization(const Mesh1D_SPtr_t& m, F func)
{
  static_assert(std::is_same<typename traits::callable_traits<F>::result_type, double>::value, "static assert error: callable type F must return a double"); 
  static_assert(traits::callable_traits<F>::num_args <= 1, "static assert error: callable type F must take <= args for Mesh1D"); 

  linops::Vector1D result(m); 

  constexpr std::size_t N = traits::callable_traits<F>::num_args; 
  if constexpr(N == 1){
    std::transform(m->cbegin(), m->cend(), result.begin(), func); 
  }
  else if constexpr(N == 0){
    result.values().setConstant( func() ); 
  }
  return result; 
}

} // end namespace linops 

#endif // Vector1D.hpp