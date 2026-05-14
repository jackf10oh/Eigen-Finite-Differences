// AutonomousCoeff.hpp 
//
// Class for coeffs of form 
// c(x,y,z) i.e. no time t given to args 
//
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_AUTONOMOUSCOEFF_H
#define FORNFDM_COEFFS_AUTONOMOUSCOEFF_H

#include "CoeffBase.hpp"

namespace fornfdm{
namespace linops{

// forward declarartion 
template<class Callable>
class AutonomousCoeff; 

namespace internal{

// Traits
template<class Callable>
struct traits_impl<fornfdm::linops::AutonomousCoeff<Callable>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = fornfdm::internal::callable_traits<Callable>::arity; 
  static constexpr bool is_timedep = false; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#include "CoeffBase.hpp"

namespace Eigen{
namespace internal{

// Traits 
template<class Callable>
struct traits<fornfdm::linops::AutonomousCoeff<Callable>> // : public traits<fornfdm::linops::CoeffBase<fornfdm::linops::AutonomousCoeff<Callable>>>
{
  typedef Eigen::DiagonalShape XprKind; 
  typedef typename Eigen::CwiseNullaryOp<fornfdm::linops::CyclicWrapper, Eigen::Matrix<fornfdm::Scalar, 1, Eigen::Dynamic>> DiagonalVectorType;
  typedef fornfdm::Scalar Scalar; 
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

}
}

namespace fornfdm{
namespace linops{  

template<class Callable>
class AutonomousCoeff : public CoeffBase<AutonomousCoeff<Callable>>
{
  private:
    // Type Defs ------------
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
    const auto& callable() const { return m_callable; }
    void setMesh(const std::shared_ptr<const fornfdm::Mesh>& m)
    {
      m_mesh_observed = m; 
      CoeffBase<AutonomousCoeff<Callable>>::setMesh_impl(m.get()); 
    }
    auto getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fornfdm::Real t){/* do nothing */}
    fornfdm::Real getTime() const { return -1.0; }
    void setTime_hooked(fornfdm::Real t){/* do nothing */}
};

} // end namespace linops 
} // end namespace fornfdm 

#endif // AutonomousCoeff.hpp 