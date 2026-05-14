// TimeDepCoeff.hpp 
//
// Class for coeffs of form 
// c(t,x,y,z) time t given to args.  
//
// JAF 4/24/2026 

#ifndef FORNFDM_COEFFS_TIMEDEPCOEFF_H
#define FORNFDM_COEFFS_TIMEDEPCOEFF_H

#include "CoeffBase.hpp"

namespace fornfdm{
namespace linops{

// forward declarartion 
template<class Callable>
class TimeDepCoeff; 

namespace internal{

// Traits
template<class Callable>
struct traits_impl<fornfdm::linops::TimeDepCoeff<Callable>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = false; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false;
  static constexpr std::size_t max_num_args_called = fornfdm::internal::callable_traits<Callable>::arity - 1; // first arg binded to time. 
  static constexpr bool is_timedep = true; 
}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#include "CoeffBase.hpp"

namespace Eigen{
namespace internal{

// Traits 
template<class Callable>
struct traits<fornfdm::linops::TimeDepCoeff<Callable>> // : public traits<fornfdm::linops::CoeffBase<fornfdm::linops::AutonomousCoeff<Callable>>>
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
class TimeDepCoeff : public CoeffBase<TimeDepCoeff<Callable>>
{
  private:
    // Type Defs ------------
    using CallableCleaned = std::remove_cv_t<std::remove_reference_t<Callable>>; 
    // Member Data ----------- 
    std::weak_ptr<const Mesh> m_mesh_observed; 
    const Mesh* m_mesh_raw; 
    typename fornfdm::internal::BindFirst<Callable> m_callable; // stores Callable + captured something that converts to fornfdm::Real. 

  public:
    // Constructor ----------
    TimeDepCoeff(Callable c)
      : CoeffBase<TimeDepCoeff<Callable>>(), m_callable(c, -1.0)
    {}
    
    // Member Functions 
    const auto& callable() const { return m_callable; }
    void setMesh(const std::shared_ptr<const Mesh>& m)
    {
      m_mesh_observed = m; 
      m_mesh_raw = m.get(); 
    }
    auto getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fornfdm::Real t)
    { 
      m_callable.captured_arg = t; 
      CoeffBase<TimeDepCoeff>::setMesh_impl(m_mesh_raw);
    }
    // fornfdm::Real getTime() const { return m_callable.captured_arg; } // don't want to do this. getTime() should always reflect last setTime that triggered updates to matrix 
    void setTime_hooked(fornfdm::Real t)
    {
      m_callable.captured_arg = t;
    }
};

} // end namespace linops 
} // end namespace fornfdm 

#endif // TimeDepCoeff.hpp 