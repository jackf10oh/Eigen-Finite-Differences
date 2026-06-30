// KroneckerEvaluator.hpp 
//
// Handles reading the CSRMatrix inside of a PartialDerivBase<> 
// into a kronecker product eigen expression  
//
// JAF 4/17/2026 

#ifndef FORNFDM_DIFFOPS_KRONECKEREVALUATOR_H
#define FORNFDM_DIFFOPS_KRONECKEREVALUATOR_H

#include "../types.hpp" // CSRMatrix
#include "../traits.hpp" // callable_traits<>

namespace fornfdm{
namespace linops{
namespace internal{

// Generic Case. Use single kronecker product to make block diagonal: I @ D 
template<class Derived, std::size_t num_args = linops::internal::traits<Derived>::max_num_args_called>
struct KroneckerEvaluator
{
  // Type Defs --------------
  typedef Derived XprType;
  typedef typename Eigen::internal::nested_eval< fornfdm::CSRMatrix, XprType::ColsAtCompileTime >::type ArgTypeNested;
  typedef typename Eigen::internal::remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 
  enum { CoeffReadCost = Eigen::internal::evaluator<fornfdm::CSRMatrix>::CoeffReadCost, Flags = Eigen::RowMajorBit };

  struct InnerIterator
  {
    // Constructors ==================================
    InnerIterator(const KroneckerEvaluator& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.m_stencil.cols() * (row_idx / eval.m_xpr.m_stencil.rows())), // no more / m_prod_before, or + row % m_prod_before. it should always be 1 in this case
      m_wrapped_it(eval.m_argImpl, row_idx % eval.m_xpr.m_stencil.rows()), // no more / m_prod_before. it should be 1 in this case. 
      m_kronImpl(eval)
    {}
    
    InnerIterator(const KroneckerEvaluator<linops::PartialDerivBase<Derived>, 0>& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.m_stencil.cols() * (row_idx / eval.m_xpr.m_stencil.rows())), // as above...  
      m_wrapped_it(eval.m_argImpl, row_idx % eval.m_xpr.m_stencil.rows()), // as above... 
      m_kronImpl(eval)
    {}

    // Member Functions ---------------------------------- 
    operator bool() const { return m_wrapped_it; }
    void operator++(){ ++m_wrapped_it; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_wrapped_it.col(); } // no more * m_prod_before. since it should be 1 in this case
    Index index() const { return m_col_offset + m_wrapped_it.col(); }
    Scalar value() const { return m_wrapped_it.value(); }

    // Member Data ----------------------------------------
    Index m_row_idx; 
    Index m_col_offset; 
    const KroneckerEvaluator& m_kronImpl; 
    typename Eigen::internal::evaluator<fornfdm::CSRMatrix>::InnerIterator m_wrapped_it; 
  }; 

  // Constructors ===================================== 
  KroneckerEvaluator(const Derived& xpr_d)
    : m_xpr(xpr_d), 
    m_rows(xpr_d.rows()), 
    m_cols(xpr_d.cols()), 
    m_argImpl(xpr_d.m_stencil)
  {}

  KroneckerEvaluator(const linops::PartialDerivBase<Derived>& xpr)
    : m_xpr(xpr.derived()), 
    m_rows(m_xpr.rows()), 
    m_cols(m_xpr.cols()), 
    m_argImpl(m_xpr.m_stencil)
  {}

  // Member Functions -------------------------------- 
  Index rows() const {return m_rows; }
  Index cols() const {return m_cols; }
  Index outerSize() const { return m_rows; }
  Index innerSize() const { return m_cols; }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Member Data ------------------------------------
  const Derived& m_xpr; 
  Index m_rows; // since .rows() / .cols() involves a multiply inside PartialDerivBase. its calculated + stored 
  Index m_cols; 
  typename Eigen::internal::evaluator<fornfdm::CSRMatrix> m_argImpl; 
};

// Specialization. Use double kronecker product into higher dimension + block diagonal: I @ D @ I
template<class Derived>
struct KroneckerEvaluator<Derived,0> 
{
  // Type Defs --------------
  typedef Derived XprType;
  typedef typename Eigen::internal::nested_eval< fornfdm::CSRMatrix, XprType::ColsAtCompileTime >::type ArgTypeNested;
  typedef typename Eigen::internal::remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 
  enum { CoeffReadCost = Eigen::internal::evaluator<fornfdm::CSRMatrix>::CoeffReadCost, Flags = Eigen::RowMajorBit };

  struct InnerIterator
  {
    // Constructors ==================================
    InnerIterator(const KroneckerEvaluator& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.m_prod_before * eval.m_xpr.m_stencil.cols() * (row_idx / (eval.m_xpr.m_prod_before * eval.m_xpr.m_stencil.rows())) + row_idx % eval.m_xpr.m_prod_before), 
      m_wrapped_it(eval.m_argImpl, (row_idx / eval.m_xpr.m_prod_before) % eval.m_xpr.m_stencil.rows()),
      m_kronImpl(eval)
    {}
    
    InnerIterator(const KroneckerEvaluator<linops::PartialDerivBase<Derived>, 0>& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.m_prod_before * eval.m_xpr.m_stencil.cols() * (row_idx / (eval.m_xpr.m_prod_before * eval.m_xpr.m_stencil.rows())) + row_idx % eval.m_xpr.m_prod_before), 
      m_wrapped_it(eval.m_argImpl, (row_idx / eval.m_xpr.m_prod_before) % eval.m_xpr.m_stencil.rows()),
      m_kronImpl(eval)
    {}

    // Member Functions ---------------------------------- 
    operator bool() const { return m_wrapped_it; }
    void operator++(){ ++m_wrapped_it; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_wrapped_it.col() * m_kronImpl.m_xpr.m_prod_before; }
    Index index() const { return m_col_offset + m_wrapped_it.col() * m_kronImpl.m_xpr.m_prod_before; }
    Scalar value() const { return m_wrapped_it.value(); }

    // Member Data ----------------------------------------
    Index m_row_idx; 
    Index m_col_offset; 
    const KroneckerEvaluator& m_kronImpl; 
    typename Eigen::internal::evaluator<fornfdm::CSRMatrix>::InnerIterator m_wrapped_it; 
  }; 

  // Constructors ===================================== 
  KroneckerEvaluator(const Derived& xpr_d)
    : m_xpr(xpr_d), 
    m_rows(xpr_d.rows()), 
    m_cols(xpr_d.cols()), 
    m_argImpl(xpr_d.m_stencil)
  {}

  KroneckerEvaluator(const linops::PartialDerivBase<Derived>& xpr)
    : m_xpr(xpr.derived()), 
    m_rows(m_xpr.rows()), 
    m_cols(m_xpr.cols()), 
    m_argImpl(m_xpr.m_stencil)
  {}

  // Member Functions -------------------------------- 
  Index rows() const {return m_rows; }
  Index cols() const {return m_cols; }
  Index outerSize() const { return m_rows; }
  Index innerSize() const { return m_cols; }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Member Data ------------------------------------
  const Derived& m_xpr; 
  Index m_rows; // since .rows() / .cols() involves a multiply inside PartialDerivBase. its calculated + stored 
  Index m_cols; 
  typename Eigen::internal::evaluator<fornfdm::CSRMatrix> m_argImpl; 
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#endif // KroneckerEvaluator.hpp