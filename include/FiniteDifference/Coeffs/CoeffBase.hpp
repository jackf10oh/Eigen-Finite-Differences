// CoeffBase.hpp 
//
// A class to represent coefficients of the form 
// c(x,y,z) or c(t,x,y,z)
// in the pde equations 
// Utt = c(x,y) * Uxx
// 
// JAF 4/24/2026 

#ifndef FDM_COEFFS_COEFFSBASE_H
#define FDM_COEFFS_COEFFSBASE_H

#include<memory>
#include<Eigen/Core> // DiagonalMatrixBase 
#include "../Diffops/Traits.hpp"
#include "../Diffops/EvaluatorBase.hpp" // Coordinate struct 

namespace fdm{
namespace linops{

// forward declaration 
template<class Derived>
class CoeffBase; 

// helper struct 
struct CyclicWrapper{
  const Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>& m_wrapped; 
  std::size_t m_mod; 
  CyclicWrapper(const Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>& v, std::size_t m)
    : m_wrapped(v), m_mod(m)
  {} 
  const fdm::Scalar& operator()(std::size_t idx) const { return m_wrapped[idx%m_mod]; }
};

namespace internal{

// traits 
template<class Derived>
struct traits_impl<fdm::linops::CoeffBase<Derived>> : public traits<Derived>{};

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm 

namespace Eigen{
namespace internal{

// Traits 
template<typename Derived>
struct traits<fdm::linops::CoeffBase<Derived>>
{
  typedef Eigen::DiagonalShape XprKind; 
  typedef typename Eigen::CwiseNullaryOp<fdm::linops::CyclicWrapper, Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>> DiagonalVectorType;
  typedef fdm::Scalar Scalar; 
  typedef typename DiagonalVectorType::StorageKind StorageKind;
  typedef typename DiagonalVectorType::StorageIndex StorageIndex;
  enum{
    RowsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
    ColsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
    MaxRowsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
    MaxColsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
    CoeffReadCost = 1, 
    Flags = Eigen::NestByRefBit
  }; 
}; 

} // end namespace internal 
} // end namespace eigen 

#include "CoeffProduct.hpp" 

namespace fdm{
namespace linops{ 

template<class Derived>
class CoeffBase: public Eigen::DiagonalBase<Derived>
{
  public:
    // Type Defs ------------------------
    typedef typename Eigen::CwiseNullaryOp<CyclicWrapper, Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>> DiagonalVectorType; 
    typedef Eigen::DiagonalShape StorageKind; 
    typedef typename DiagonalVectorType::StorageIndex StorageIndex; 
    using Base = Eigen::DiagonalBase<Derived>; 

  private:
    // Member Data ------------------------ 
    typename Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic> m_diagonal; 
    StorageIndex m_prod_after; 
    DiagonalVectorType m_cyclic_wrapper; 

  public:
    // Constructors 
    CoeffBase()
      : m_diagonal(0), 
      m_prod_after(0), 
      m_cyclic_wrapper(1,0,CyclicWrapper(m_diagonal, m_prod_after))
    {} 

    // Member Functions -------------------
    using Base::derived; 
    // using Base::const_derived; ????  
    const auto& diagonal() const { return m_cyclic_wrapper; }
    const auto& callable() const { return derived().callable(); }

    // Operators ------------
    template<class RhsDeriv, typename = std::enable_if_t<linops::internal::is_partialderiv_crtp<RhsDeriv>::value>>
    auto operator*(RhsDeriv&& rhs) & 
    {
      return CoeffProduct<Derived&, RhsDeriv>(derived(), std::forward<RhsDeriv>(rhs)); 
    } 

    template<class RhsDeriv, typename = std::enable_if_t<linops::internal::is_partialderiv_crtp<RhsDeriv>::value>>
    auto operator*(RhsDeriv&& rhs) && 
    {
      return CoeffProduct<Derived, RhsDeriv>(std::move(derived()), std::forward<RhsDeriv>(rhs)); 
    }

  public: // TODO encapsulate this 
    // Implementations ---------------------
    void setMesh_impl(const Mesh* m)
    {
      using traits_t = fdm::linops::internal::traits<Derived>; 
      m_prod_after = m->sizesMiddleProduct(traits_t::max_num_args_called+1, m->numDims());
      
      std::size_t end = m->sizesMiddleProduct(0, traits_t::max_num_args_called+1);
      m_diagonal.resize(end); 
      for(std::size_t idx=0; idx<end; ++idx)
      {
        fdm::Coordinate<traits_t::max_num_args_called> coord(m,idx);
        m_diagonal[idx] = coord.apply(callable());  
      }
      // placement new shenanigans
      new (&m_cyclic_wrapper) DiagonalVectorType(1,end*m_prod_after,CyclicWrapper(m_diagonal, end)); 
    }
};

} // end namespace linops
} // end namespace fdm 

#endif // CoeffBase.hpp 