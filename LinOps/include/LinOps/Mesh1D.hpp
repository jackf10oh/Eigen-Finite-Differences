// Mesh1D.hpp
//
// 1D Mesh classes for rectangular domains
// representing Domain x0 , ..., xN 
// that Function Discretizations / Linear Operators (integration/derivative) work on
//
// JAF 12/5/2025

#ifndef MESH1D_H
#define MESH1D_H

#include<cstdint>
#include<vector>
#include<memory> 
#include<type_traits>

namespace linops{

// forward declaration -> aliases
class Mesh1D; 
using SharedConstMesh1D = std::shared_ptr<const Mesh1D>;
using WeakConstMesh1D = std::weak_ptr<const Mesh1D>;

class Mesh1D
{
  protected:
    // Member Data ------------------------------------------------------------------------
    std::vector<double> m_vals;
  public:
    // Constructors / Destructor =========================================================================
    // uniformly on interval [x1,x2] 
    Mesh1D(double x1=0.0, double x2=1.0, std::size_t n_steps=11)
      : m_vals(n_steps)
    {
      if(n_steps < 2) throw std::invalid_argument("must have # of steps >= 2 in mesh");
      if(x1 >= x2) throw std::invalid_argument("must have x1 < x2 in mesh");
      double dx = (x2 - x1) / (n_steps-1); 
      for(std::size_t i=0; i<n_steps-1; i++) m_vals[i] = x1 + i*dx; 
      m_vals[n_steps-1] = x2; 
    }
    // from std::vector
    Mesh1D(const std::vector<double>& vals_init) : m_vals(vals_init){}; 
    Mesh1D(std::vector<double>&& vals_init) : m_vals(std::move(vals_init)){}; 
    // copy constructors
    Mesh1D(const Mesh1D& other): m_vals(other.m_vals){}; 

    // destructors 
    virtual ~Mesh1D()=default; 

    // Member Functions ========================================================================
    // size of mesh ---------------------------
    std::size_t size() const { return m_vals.size(); }
    
    // indexing ----------------------------
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ return m_vals.at(i); }
    const double& at(std::size_t i) const { return m_vals.at(i); }
    
    // forward iterators --------------------------
    auto begin(){ return m_vals.begin(); }
    auto end(){ return m_vals.end(); }
    
    // forward citerators ------------------------
    auto cbegin() const { return m_vals.cbegin(); }
    auto cend() const { return m_vals.cend(); }
    
    // reverse iterators --------------------
    auto rbegin(){ return m_vals.rbegin(); }
    auto rend(){ return m_vals.rend(); }

    // reverse citerators -------------------------------------------
    auto crbegin() const { return m_vals.crbegin(); }
    auto crend() const { return m_vals.crend(); }  
};

// quick helper to wrap std::make_shared<Mesh1D> or other Mesh1D derived types  
template<typename... Args> 
auto make_Mesh1D(Args... args)
{
  return std::make_shared<const Mesh1D>(args...); 
}

} // end namespace linops 

#endif // Mesh1D.hpp