// CoeffBase.hpp 
//
// A class to represent coefficients of the form 
// c(x,y,z) or c(t,x,y,z)
// in the pde equations 
// Utt = c(x,y) * Uxx
// 
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_COEFFSBASE_H
#define FORNFDM_COEFFS_COEFFSBASE_H

#include<memory>
#include<Eigen/Core> // DiagonalMatrixBase 
#include "../Diffops/Traits.hpp"
#include "../Diffops/EvaluatorBase.hpp" // Coordinate struct 
#include "../Coordinate.hpp"

namespace fornfdm{
namespace linops{

// forward declaration 
template<class Derived>
class CoeffBase; 

// helper struct 
struct CyclicWrapper{
  const Eigen::Matrix<fornfdm::Scalar, 1, Eigen::Dynamic>& m_wrapped; 
  std::size_t m_mod; 
  CyclicWrapper(const Eigen::Matrix<fornfdm::Scalar, 1, Eigen::Dynamic>& v, std::size_t m)
    : m_wrapped(v), m_mod(m)
  {} 
  const fornfdm::Scalar& operator()(std::size_t idx) const { return m_wrapped[idx%m_mod]; }
};

namespace internal{

// traits 
template<class Derived>
struct traits_impl<fornfdm::linops::CoeffBase<Derived>> : public traits<Derived>{};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

namespace Eigen{
namespace internal{

// Traits ------ 
template<typename Derived>
struct traits<fornfdm::linops::CoeffBase<Derived>> : public traits<Derived>
{
  enum {Flags = traits<Derived>::Flags}; 
}; 

} // end namespace internal 
} // end namespace Eigen 

#include "CoeffProduct.hpp" 

namespace fornfdm{
namespace linops{ 

template<class Derived>
class CoeffBase: public Eigen::DiagonalBase<Derived>
{
  public:
    // Type Defs ------------------------
    typedef typename Eigen::DiagonalBase<Derived> Base; 
    typedef typename Eigen::CwiseNullaryOp<CyclicWrapper, Eigen::Matrix<fornfdm::Scalar, 1, Eigen::Dynamic>> DiagonalVectorType; 
    // typedef typename Eigen::internal::traits<Derived>::DiagonalVectorType DiagonalVectorType;
    typedef typename DiagonalVectorType::Scalar Scalar;
    typedef typename DiagonalVectorType::RealScalar RealScalar;
    typedef typename Eigen::internal::traits<Derived>::StorageKind StorageKind;
    typedef typename Eigen::internal::traits<Derived>::StorageIndex StorageIndex;
    enum {
      RowsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
      ColsAtCompileTime = DiagonalVectorType::SizeAtCompileTime,
      MaxRowsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
      MaxColsAtCompileTime = DiagonalVectorType::MaxSizeAtCompileTime,
      IsVectorAtCompileTime = 0,
      Flags = Eigen::NoPreferredStorageOrderBit
    };
    typedef Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime, 0, MaxRowsAtCompileTime, MaxColsAtCompileTime> DenseMatrixType;
    typedef DenseMatrixType DenseType;
    typedef Eigen::DiagonalMatrix<Scalar,DiagonalVectorType::SizeAtCompileTime,DiagonalVectorType::MaxSizeAtCompileTime> PlainObject;

  protected:
    // Member Data ------------------------ 
    typename Eigen::Matrix<fornfdm::Scalar, 1, Eigen::Dynamic> m_diagonal; 
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
    void setMesh(const std::shared_ptr<const fornfdm::Mesh>& m){ return derived().setMesh(m); }
    auto getMesh() const { return derived().getMesh(); }
    void setTime(fornfdm::Real t){ derived().setTime(t); }
    fornfdm::Real getTime() const { return derived().getTime(); }
    const auto& toEigen() const { return *static_cast<const Eigen::DiagonalBase<Derived>*>(this); }
    using Base::derived; 
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
};

} // end namespace linops
} // end namespace fornfdm 

#endif // CoeffBase.hpp 