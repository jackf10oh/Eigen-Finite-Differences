// AutonomousCoeff.hpp 
//
// Class for coeffs of form 
// c(x,y,z) i.e. no time t given to args 
//
// JAF 4/24/2026 

#ifndef FDM_COEFFS_AUTONOMOUSCOEFF_H
#define FDM_COEFFS_AUTONOMOUSCOEFF_H

#include "CoeffBase.hpp"

namespace fdm{
namespace linops{

// forward declarartion 
template<class Callable>
class AutonomousCoeff; 

namespace internal{

// Traits
template<class Callable>
struct traits_impl<fdm::linops::AutonomousCoeff<Callable>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_num_args_called = fdm::internal::callable_traits<Callable>::arity; 
  static constexpr bool is_timedep = false; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm 

#include "CoeffBase.hpp"

namespace Eigen{
namespace internal{

// Traits 
template<class Callable>
struct traits<fdm::linops::AutonomousCoeff<Callable>> // : public traits<fdm::linops::CoeffBase<fdm::linops::AutonomousCoeff<Callable>>>
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

}
}

namespace fdm{
namespace linops{  

template<class Callable>
class AutonomousCoeff : public CoeffBase<AutonomousCoeff<Callable>>
{
  private:
    // Type Defs ------------
    using CallableCleaned = std::remove_cv_t<std::remove_reference_t<Callable>>; 
    // Member Data ----------- 
    std::weak_ptr<const fdm::Mesh> m_mesh_observed; 
    CallableCleaned m_callable; 

  public:
    // Constructor ----------
    AutonomousCoeff(Callable c)
      : CoeffBase<AutonomousCoeff<Callable>>(), m_callable(c)
    {}
    
    // Member Functions 
    const auto& callable() const { return m_callable; }
    void setMesh(const std::shared_ptr<const fdm::Mesh>& m)
    {
      m_mesh_observed = m; 
      CoeffBase<AutonomousCoeff<Callable>>::setMesh_impl(m.get()); 
    }
    auto getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fdm::Real t){/* do nothing */}
    fdm::Real getTime() const { return -1.0; }
    void setTime_hooked(fdm::Real t){/* do nothing */}
};

} // end namespace linops 
} // end namespace fdm 

#endif // AutonomousCoeff.hpp 