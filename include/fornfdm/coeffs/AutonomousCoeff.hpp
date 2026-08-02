// AutonomousCoeff.hpp 
//
// Class for coeffs of form 
// c(x,y,z) i.e. no time t given to args 
//
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_AUTONOMOUSCOEFF_H
#define FORNFDM_COEFFS_AUTONOMOUSCOEFF_H

#include<cstdint>
#include<Eigen/Core> // DiagonalMatrixBase 
#include "../types.hpp"
#include "../traits.hpp"
#include "../diffops/traits.hpp"
#include "../diffops/functors.hpp"
#include "../Mesh.hpp"
#include "../Coordinate.hpp"
#include "CoeffBase.hpp"

namespace fornfdm{
namespace linops{

// forward declarartion 
template<class Callable>
class AutonomousCoeff; 

namespace internal{

// traits
template<class Callable>
struct traits_impl<fornfdm::linops::AutonomousCoeff<Callable>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = fornfdm::internal::callable_traits<Callable>::arity; 
  static constexpr bool is_timedep = false; 
  using orders = std::index_sequence<>;
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

namespace Eigen{
namespace internal{

// traits 
template<class Callable>
struct traits<fornfdm::linops::AutonomousCoeff<Callable>> 
  : public traits<Eigen::CwiseNullaryOp<fornfdm::linops::CyclicWrapper, fornfdm::Vector>>
{
  typedef typename Eigen::CwiseNullaryOp<fornfdm::linops::CyclicWrapper, fornfdm::Vector> DiagonalVectorType;
  typedef DiagonalShape StorageKind;
  enum {
    Flags = LvalueBit | NoPreferredStorageOrderBit
  };
};

} // end namespace internal
} // end namespace Eigen

namespace fornfdm{
namespace linops{  

template<class Callable>
class AutonomousCoeff : public CoeffBase<AutonomousCoeff<Callable>>
{
  public:
    // Type Defs ------------
    typedef typename Eigen::internal::traits<AutonomousCoeff>::DiagonalVectorType DiagonalVectorType;
    typedef const AutonomousCoeff& Nested;
    typedef fornfdm::Scalar Scalar;
    typedef typename Eigen::internal::traits<AutonomousCoeff>::StorageKind StorageKind;
    typedef typename Eigen::internal::traits<AutonomousCoeff>::StorageIndex StorageIndex;
  private:
    using CallableCleaned = std::remove_cv_t<std::remove_reference_t<Callable>>; 
    // Member Data ----------- 
    std::weak_ptr<const fornfdm::Mesh> m_mesh_observed; 
    CallableCleaned m_callable; 

  public:
    // Constructor ----------
    AutonomousCoeff(Callable c)
      : CoeffBase<AutonomousCoeff<Callable>>(), m_callable(c)
    {}
    
    // Member Functions 
    using CoeffBase<AutonomousCoeff<Callable>>::diagonal; 
    const auto& callable() const { return m_callable; }
    void setMesh(const std::shared_ptr<const fornfdm::Mesh>& m)
    {
      m_mesh_observed = m; 
      using traits_t = fornfdm::linops::internal::traits<AutonomousCoeff<Callable>>; 
      this->m_prod_after = m->sizesMiddleProduct(traits_t::max_arity, m->numDims());
      
      std::size_t end = m->sizesMiddleProduct(0, traits_t::max_arity);
      this->m_diagonal.resize(end); 
      for(std::size_t idx=0; idx<end; ++idx)
      {
        fornfdm::Coordinate<traits_t::max_arity> coord(m.get(),idx);
        this->m_diagonal[idx] = coord.apply(m_callable);  
      }
      // placement new shenanigans
      new (&(this->m_cyclic_wrapper)) typename CoeffBase<AutonomousCoeff<Callable>>::DiagonalVectorType(end*(this->m_prod_after),1,CyclicWrapper(this->m_diagonal, end)); 
    }
    auto getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fornfdm::Real t){/* do nothing */}
    fornfdm::Real getTime() const { return -1.0; }
};

} // end namespace linops 
} // end namespace fornfdm 

#endif // AutonomousCoeff.hpp 