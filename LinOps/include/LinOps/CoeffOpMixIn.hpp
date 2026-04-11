// CoeffOpMixIn.hpp
//
// going to need some way to model equations like 
// Ut = c(t,x) * Uxx + c(t,x) * Ux
//
// CoeffOpMixIn.hpp

#ifndef COEFFOP_H
#define COEFFOP_H 

#include<Eigen/Core> 
#include "LinOpTraits.hpp" 
#include "LinOpMixIn.hpp"
#include "LinOpBase.hpp"

namespace linops{

template<typename Derived, typename Callable>
class CoeffOpMixIn : public LinOpMixIn<Derived>
{

  public:
    // Type Defs --------------------------------
    using DerivedT = Derived; // So LinOpMixIn case access Derived type of grandchild classes
    struct is_coeff_flag{}; // flag type for any derived class to be picked up by is_coeffop_crtp<>::value trait  
    static constexpr std::size_t numArgs = traits::callable_traits<Callable>::num_args; 
  protected:
    // Member Data -------------------------------------- 
    Callable m_callable; 
    Eigen::VectorXd m_diag; 
    std::size_t m_prod_after; 

  public:
    // Constructors + Destructor ========================================================= 
    // default 
    CoeffOpMixIn()=delete;
    
    // from callable function f 
    CoeffOpMixIn(Callable f) 
      : m_callable(std::move(f)), m_diag(), m_prod_after() 
    {}

    // copy 
    CoeffOpMixIn(const CoeffOpMixIn& other)=default; 

    // destructor 
    ~CoeffOpMixIn()=default; 

    // Member Functions ============================================================
    auto asMatrix() const 
    {
      constexpr std::size_t n = linops::traits::callable_traits<Callable>::num_args; 
      if constexpr(n==0)
      {
        auto repeated = Eigen::VectorXd::NullaryExpr( m_prod_after, [val = m_diag[0]](std::size_t idx){ return val; });  
        return repeated.asDiagonal(); 
      }
      else
      {
        auto repeated = Eigen::VectorXd::NullaryExpr(m_diag.size() * m_prod_after, [&](std::size_t idx){ return m_diag[idx % m_diag.size()]; }); 
        return  repeated.asDiagonal(); 
      }
    }

    // Must implement! 
    void setMesh1D_impl(const SharedConstMesh1D& m)
    {
      static_cast<DerivedT*>(this)->setMesh1D_impl(m); 
    }
  
    // Must implement! 
    void setMeshXD_impl(const SharedConstMeshXD& m)
    {
      static_cast<DerivedT*>(this)->setMeshXD_impl(m); 
    }
    
    const Callable& callable() const { return m_callable; }

    void fillDiagonal(const linops::SharedConstMesh1D& m)
    {
      constexpr std::size_t n = linops::traits::callable_traits<Callable>::num_args; 
      if constexpr(n==0){
        m_prod_after = m->size(); 
        m_diag.resize(1); 
        m_diag[0] = m_callable(); 
      }
      else if constexpr(n==1){
        m_prod_after = 1;
        std::size_t s = m->size();  
        m_diag.resize(s); 
        for(auto i=0; i<s; ++i) m_diag[i] = m_callable((*m)[i]); 
      }
      else{
        throw std::runtime_error("CoeffOpMixin error. fillDiagonal called on Mesh1D when callable had > 1 arguments"); 
      }
    }

    void fillDiagonal(const linops::SharedConstMeshXD& m)
    {
      constexpr std::size_t n = linops::traits::callable_traits<Callable>::num_args; 
      if(n > m->numDims()) throw std::runtime_error("CoeffOpMixin error. fillDiagonal called on MeshXD when callale had num args > numDims in mesh");

      std::size_t end = m->sizesMiddleProduct(0,n); 
      m_diag.resize(end); 
      m_prod_after = m->sizesMiddleProduct(n, m->numDims());
      
      std::array<double,n> coords; 
      std::array<std::size_t,n> prods_arr;
      std::size_t rolling = 1; // rolling product of Mesh1D sizes... 
      for(std::size_t ith_dim=0; ith_dim < n; ++ith_dim)
      {
        rolling *= m->sizeOfDim(ith_dim); 
         prods_arr[ith_dim] = rolling; 
      }
      
      for(std::size_t flat_idx=0; flat_idx<end; ++flat_idx)
      {
        for(std::size_t ith_dim=0; ith_dim < n; ++ith_dim)
        {
          coords[ith_dim] = (*(m->getMesh1D(ith_dim)))[flat_idx % prods_arr[ith_dim]]; 
        }
        m_diag[flat_idx] = std::apply(m_callable, coords); 
      }
    }

    // Operators -----------------------------------------------------------------
    // composition c * L  produces composition
    template<typename RHS>
    auto operator*(RHS&& rhs)
    {
      static_assert(!traits::is_coeffop_crtp<RHS>::value,"Coefficients are meant to multiply c*L for L linear operator. not another Coefficient. a*b*L should be written as 1 functions");
      return LinOpMixIn<Derived>::compose(std::forward<RHS>(rhs));
    };
    
    // delting a ton of operators out of LinOpMixIn 
    // so that you can't make expressions of coefficients .......
    auto operator-()=delete;
    template<typename RHS>
    auto operator-(RHS&& rhs)=delete;
    template<typename RHS>
    auto operator+(RHS&& rhs)=delete; 
    template<typename RHS>
    auto compose(RHS&& InnerOp)=delete; 
    auto left_scalar_mult_impl(double c)=delete; 
};

} // end namespace linops 

#endif // CoeffOpMixIn.hpp 