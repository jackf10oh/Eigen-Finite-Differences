// TimeDepCoeff.hpp 
//
// Class for coeffs of form 
// c(t,x,y,z) time t given to args.  
//
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_TIMEDEPCOEFF_H
#define FORNFDM_COEFFS_TIMEDEPCOEFF_H

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
class TimeDepCoeff; 

namespace internal{

// traits
template<class Callable>
struct traits_impl<fornfdm::linops::TimeDepCoeff<Callable>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = fornfdm::internal::callable_traits<Callable>::arity - 1; // first argument is time
  static constexpr bool is_timedep = true; 
  using orders = std::index_sequence<>;
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

namespace Eigen{
namespace internal{

// traits 
template<class Callable>
struct traits<fornfdm::linops::TimeDepCoeff<Callable>> 
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
class TimeDepCoeff : public CoeffBase<TimeDepCoeff<Callable>>
{
  public:
    // Type Defs ------------
    typedef typename Eigen::internal::traits<TimeDepCoeff>::DiagonalVectorType DiagonalVectorType;
    typedef const TimeDepCoeff& Nested;
    typedef fornfdm::Scalar Scalar;
    typedef typename Eigen::internal::traits<TimeDepCoeff>::StorageKind StorageKind;
    typedef typename Eigen::internal::traits<TimeDepCoeff>::StorageIndex StorageIndex;
  private:
    using CallableCleaned = std::remove_cv_t<std::remove_reference_t<Callable>>; 
    // Member Data ----------- 
    std::weak_ptr<const Mesh> m_mesh_observed; 
    const Mesh* m_mesh_raw;
    fornfdm::Real m_current_time;  
    CallableCleaned m_callable;

  public:
    // Constructor ----------
    TimeDepCoeff(Callable c)
      : CoeffBase<TimeDepCoeff<Callable>>(), m_callable(c)
    {}
    
    // Member Functions 
    using CoeffBase<TimeDepCoeff<Callable>>::diagonal; 
    const auto& callable() const { return m_callable; }
    void setMesh(const std::shared_ptr<const Mesh>& m)
    {
      m_mesh_observed = m; 
      m_mesh_raw = m.get(); 
    }
    auto getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fornfdm::Real t)
    { 
      using traits_t = fornfdm::linops::internal::traits<TimeDepCoeff<Callable>>; 
      this->m_prod_after = m_mesh_raw->sizesMiddleProduct(traits_t::max_arity, m_mesh_raw->numDims());
      
      std::size_t end = m_mesh_raw->sizesMiddleProduct(0, traits_t::max_arity);
      this->m_diagonal.resize(end); 
      for(std::size_t idx=0; idx<end; ++idx)
      {
        fornfdm::Coordinate<traits_t::max_arity> coord(m_mesh_raw,idx);
        this->m_diagonal[idx] = coord.applyBindFirst(m_callable, t);  
      }
      // placement new shenanigans
      new (&(this->m_cyclic_wrapper)) typename CoeffBase<TimeDepCoeff<Callable>>::DiagonalVectorType(end*(this->m_prod_after),1,CyclicWrapper(this->m_diagonal, end)); 
    }
    fornfdm::Real getTime() const { return m_current_time; }
};

} // end namespace linops 
} // end namespace fornfdm 

#endif // TimeDepCoeff.hpp 