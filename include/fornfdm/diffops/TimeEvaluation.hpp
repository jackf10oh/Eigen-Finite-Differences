// TimeEvaluation.hpp
//
// Eigen Expression for the single shot evaluation
// of a diffops expression at T=t0 without updating internal state. 
//
// JAF 7/19/2026

#ifndef FORNFDM_DIFFOPS_TIMEEVALUATION_H
#define FORNFDM_DIFFOPS_TIMEEVALUATION_H

#include "../Mesh.hpp"
#include "traits.hpp"
#include "EvaluatorBase.hpp"
#include "EigenEvaluatorImpl.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"

namespace fornfdm{
namespace linops{

// Forward Declaration ----------
template<class ArgType, class>
class TimeEvaluation;

namespace internal{

// diffops internal traits ------------- 
template<class ArgType>
struct traits_impl<TimeEvaluation<ArgType>> : traits<ArgType>{};

} // end namespace internal 

template<class ArgType>
struct TimeEvaluation<ArgType> : public Eigen::SparseMatrixBase<TimeEvaluation<ArgType>>
{
  // Type Defs ----------------
  typedef typename fornfdm::linops::TimeEvaluation<ArgType> Nested;

  // Member Data ------------------ 
  const ArgType& m_arg;
  const fornfdm::Real m_time; 

  // Constructors ----------------- 
  TimeEvaluation(const ArgType& arg, fornfdm::Real t)
    : m_arg(arg), m_time(t)
  {}

  // Member Functions --------------- 
  auto rows() const { return m_arg.rows(); } 
  auto cols() const { return m_arg.cols(); } 
  auto nonZerosEstimate() const { return m_arg.nonZerosEstimate(); }
};

} // end namespace linops
} // end namespace fornfdm

namespace Eigen{
namespace internal{

// Eigen's internal traits ----------
template<class ArgType>
struct traits<fornfdm::linops::TimeEvaluation<ArgType>>
{
  typedef fornfdm::Scalar Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Sparse StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = Eigen::RowMajorBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<class ArgType>
struct evaluator<fornfdm::linops::TimeEvaluation<ArgType>> 
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::TimeEvaluation<ArgType>>
{
  using XprType = fornfdm::linops::TimeEvaluation<ArgType>; 
  using Impl = typename fornfdm::linops::internal::EigenEvaluatorImpl<XprType>; 
  using InnerIterator = typename Impl::InnerIterator; 
  evaluator(const XprType& xpr)
    : Impl(xpr)
  {}
}; 

  } // end namespace internal
} // end namespace Eigen

#endif // TimeEvaluation.hpp

