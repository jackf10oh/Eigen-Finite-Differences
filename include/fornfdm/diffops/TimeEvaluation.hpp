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
#include "PartialDerivBase.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"

namespace fornfdm{
  namespace linops{

// Forward Declaration ----------
template<class ArgType, std::size_t max_args> // also in PartialDerivBase.hpp with defaults for max_args
class TimeEvaluation;

namespace internal{

// diffops internal traits ------------- 
template<class ArgType>
struct traits_impl<TimeEvaluation<ArgType>> : traits<ArgType>{};

} // end namespace internal 

  } // end namespace linops
} // end namespace fornfdm

namespace Eigen{
  namespace internal{

// Eigen's internal traits ----------
template<class ArgType, std::size_t max_args>
struct traits<fornfdm::linops::TimeEvaluation<ArgType, max_args>>
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

  } // end namespace internal  
} // end namespace Eigen

namespace fornfdm{
  namespace linops{

// The definition -----------
template<class ArgType, std::size_t max_args>
class TimeEvaluation : public Eigen::SparseMatrixBase<TimeEvaluation<ArgType,max_args>>
{
  public:
    // Friends ------------------ 
    friend Eigen::internal::evaluator<TimeEvaluation<ArgType>>;

    // Type Defs ----------------
    typedef typename fornfdm::linops::TimeEvaluation<ArgType,max_args> Nested;

  private:
    // Member Data ------------------ 
    const ArgType& m_arg;
    const Mesh* m_mesh;
    const fornfdm::Real m_time; 
    std::size_t m_size;

  public:
    // Constructors ----------------- 
    TimeEvaluation(const ArgType& arg, const fornfdm::Mesh* mesh, fornfdm::Real t)
      : m_arg(arg), m_mesh(mesh), m_time(t)
    { m_size = m_mesh->sizesProduct(); }

    TimeEvaluation(const linops::PartialDerivBase<ArgType>& pde_base, const fornfdm::Mesh* mesh, fornfdm::Real t)
      : m_arg(pde_base), m_mesh(mesh), m_time(t)
    { m_size = m_mesh->sizesProduct();}

    // Member Functions --------------- 
    std::size_t rows() const { return m_size; } 
    std::size_t cols() const { return m_size; } 
};

  } // end namespace linops
} // end namespace fornfdm

namespace Eigen{
  namespace internal{

template<class ArgType, std::size_t max_args>
struct evaluator<fornfdm::linops::TimeEvaluation<ArgType,max_args>>
{
  // Type Defs 
  typedef typename traits<fornfdm::linops::TimeEvaluation<ArgType,max_args>>::StorageIndex StorageIndex; 
  typedef typename traits<fornfdm::linops::TimeEvaluation<ArgType,max_args>>::Scalar Scalar; 
  typedef typename fornfdm::linops::TimeEvaluation<ArgType,max_args> Nested;

  // Flags -------
  enum { CoeffReadCost = evaluator<ArgType>::CoeffReadCost, Flags = Eigen::RowMajor };

  // Member Data ---------------- 
  using traits_t = fornfdm::linops::internal::traits<ArgType>;
  Nested m_xpr;
  fornfdm::linops::internal::Evaluator<ArgType> m_eval;
  StorageIndex m_prod_before;
  StorageIndex m_stencil_size;
  StorageIndex m_prod_after;

  struct InnerIterator
  {
    // Member Data -----------
    const evaluator& m_eval; 
    std::size_t m_idx;
    std::size_t m_row_idx;
    std::size_t m_offset; 
    std::size_t m_stride; 
    typename fornfdm::linops::internal::Evaluator<ArgType>::Row m_row;

    // Constructor --------
    InnerIterator(const evaluator& eval, std::size_t row_index)
      : m_eval(eval),
      m_idx(0), 
      m_row_idx(row_index),
      m_offset( m_eval.m_prod_before * m_eval.m_stencil_size * (m_row_idx / (m_eval.m_prod_before * m_eval.m_stencil_size)) + (m_row_idx % m_eval.m_prod_before) ),
      m_row(eval.m_eval, eval.m_xpr.m_mesh, (row_index/eval.m_prod_before)%(eval.m_stencil_size), row_index, eval.m_xpr.m_time)
    {}
    
    // Member Functions --------------
    operator bool() const { return (m_idx != m_row.size()); }
    void operator++(){ ++m_idx; }
    StorageIndex row() const { return m_row_idx; }
    StorageIndex col() const { return m_offset + m_row.index(m_idx) * m_eval.m_prod_before; }
    StorageIndex index() const { return m_offset + m_row.index(m_idx) * m_eval.m_prod_before; }
    Scalar value() const { return m_row.value(m_idx); }
  };

  // Constructor ---------------- 
  evaluator(const fornfdm::linops::TimeEvaluation<ArgType,max_args>& xpr)
    : m_xpr(xpr), m_eval(xpr.m_arg)
  {
    m_prod_before = m_xpr.m_mesh->sizesMiddleProduct(0,traits_t::direction);
    m_prod_after = m_xpr.m_mesh->sizesMiddleProduct(traits_t::direction+1,m_xpr.m_mesh->numDims());
    m_stencil_size = m_xpr.m_mesh->sizeOfDim(traits_t::direction);
  }

  // member Functions ------------------ 
  std::size_t rows() const { return m_xpr.rows(); }
  std::size_t cols() const { return m_xpr.cols(); }
  Index outerSize() const { return  m_xpr.rows(); }
  Index innerSize() const { return  m_xpr.cols(); }
  Index nonZerosEstimate() const { return m_prod_before * m_eval.totalNonZeros(m_xpr.m_mesh->getAxis(traits_t::direction)) * m_prod_after; }
};

  } // end namespace internal
} // end namespace Eigen

#endif // TimeEvaluation.hpp

