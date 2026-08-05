// EigenEvaluatorImpl.hpp 
//
// Handles reading the CSRMatrix inside of a PartialDerivBase<> 
// into a kronecker product eigen expression
//
// JAF 4/17/2026 

#ifndef FORNFDM_DIFFOPS_EIGENEVALUATORIMPL_H
#define FORNFDM_DIFFOPS_EIGENEVALUATORIMPL_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include "../types.hpp" // CSRMatrix
#include "../traits.hpp" // callable_traits<>
#include "PartialDerivBase.hpp"

namespace fornfdm{
namespace linops{
namespace internal{

// empty class to specialize
template<class Derived, class = void>
struct EigenEvaluatorImpl{};

// ================================================================== 
// LeftKroneckerTag + TimeDepLeftKroneckerTag. 
// Use single kronecker product to make block diagonal: I @ D
// ==================================================================

template<class Derived>
struct EigenEvaluatorImpl<Derived, std::enable_if_t<std::is_base_of_v<LeftKroneckerTag, typename map_to_base_tag<Derived>::type>>>
  : public Eigen::internal::evaluator_base<Derived>
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
    InnerIterator(const EigenEvaluatorImpl& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.getStencil().cols() * (row_idx / eval.m_xpr.getStencil().rows())),
      m_wrapped_it(eval.m_argImpl, row_idx % eval.m_xpr.getStencil().rows()),  
      m_kronImpl(eval)
    {}
    
    InnerIterator(const EigenEvaluatorImpl<linops::PartialDerivBase<Derived>>& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.getStencil().cols() * (row_idx / eval.m_xpr.getStencil().rows())), // as above...  
      m_wrapped_it(eval.m_argImpl, row_idx % eval.m_xpr.getStencil().rows()), // as above... 
      m_kronImpl(eval)
    {}

    // Member Functions ---------------------------------- 
    operator bool() const { return m_wrapped_it; }
    void operator++(){ ++m_wrapped_it; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_wrapped_it.col(); } 
    Index index() const { return m_col_offset + m_wrapped_it.col(); }
    Scalar value() const { return m_wrapped_it.value(); }

    // Member Data ----------------------------------------
    Index m_row_idx; 
    Index m_col_offset; 
    const EigenEvaluatorImpl& m_kronImpl; 
    typename Eigen::internal::evaluator<fornfdm::CSRMatrix>::InnerIterator m_wrapped_it; 
  }; 

  // Constructors ===================================== 
  EigenEvaluatorImpl(const Derived& xpr_d)
    : m_xpr(xpr_d), 
    m_rows(xpr_d.rows()), 
    m_cols(xpr_d.cols()), 
    m_argImpl(xpr_d.getStencil())
  {}

  EigenEvaluatorImpl(const linops::PartialDerivBase<Derived>& xpr)
    : m_xpr(xpr.derived()), 
    m_rows(m_xpr.rows()), 
    m_cols(m_xpr.cols()), 
    m_argImpl(m_xpr.getStencil())
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

// ==================================================================
// DoubleKroneckerTag + TimeDepDoubleKroneckerTag
// Use double kronecker product into higher dimension + block diagonal: I @ D @ I
// ==================================================================

template<class Derived>
struct EigenEvaluatorImpl<Derived, std::enable_if_t<std::is_base_of_v<DoubleKroneckerTag, typename map_to_base_tag<Derived>::type>>>
  : public Eigen::internal::evaluator_base<Derived>
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
    InnerIterator(const EigenEvaluatorImpl& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.getProductBefore() * eval.m_xpr.getStencil().cols() * (row_idx / (eval.m_xpr.getProductBefore() * eval.m_xpr.getStencil().rows())) + row_idx % eval.m_xpr.getProductBefore()), 
      m_wrapped_it(eval.m_argImpl, (row_idx / eval.m_xpr.getProductBefore()) % eval.m_xpr.getStencil().rows()),
      m_kronImpl(eval)
    {}
    
    InnerIterator(const EigenEvaluatorImpl<linops::PartialDerivBase<Derived>>& eval, Index row_idx)
      : m_row_idx(row_idx), 
      m_col_offset(eval.m_xpr.getProductBefore() * eval.m_xpr.getStencil().cols() * (row_idx / (eval.m_xpr.getProductBefore() * eval.m_xpr.getStencil().rows())) + row_idx % eval.m_xpr.getProductBefore()), 
      m_wrapped_it(eval.m_argImpl, (row_idx / eval.m_xpr.getProductBefore()) % eval.m_xpr.getStencil().rows()),
      m_kronImpl(eval)
    {}

    // Member Functions ---------------------------------- 
    operator bool() const { return m_wrapped_it; }
    void operator++(){ ++m_wrapped_it; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_wrapped_it.col() * m_kronImpl.m_xpr.getProductBefore(); }
    Index index() const { return m_col_offset + m_wrapped_it.col() * m_kronImpl.m_xpr.getProductBefore(); }
    Scalar value() const { return m_wrapped_it.value(); }

    // Member Data ----------------------------------------
    Index m_row_idx; 
    Index m_col_offset; 
    const EigenEvaluatorImpl& m_kronImpl; 
    typename Eigen::internal::evaluator<fornfdm::CSRMatrix>::InnerIterator m_wrapped_it; 
  }; 

  // Constructors ===================================== 
  EigenEvaluatorImpl(const Derived& xpr_d)
    : m_xpr(xpr_d), 
    m_rows(xpr_d.rows()), 
    m_cols(xpr_d.cols()), 
    m_argImpl(xpr_d.getStencil())
  {}

  EigenEvaluatorImpl(const linops::PartialDerivBase<Derived>& xpr)
    : m_xpr(xpr.derived()), 
    m_rows(m_xpr.rows()), 
    m_cols(m_xpr.cols()), 
    m_argImpl(m_xpr.getStencil())
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

// ==================================================================
// StoredWeightsTag + TimeDepStoredWeightsTag
// ==================================================================

// with direction == 0
template<class Derived>
struct EigenEvaluatorImpl<
  Derived, 
  std::enable_if_t<
    (std::is_base_of_v<StoredWeightsTag, typename map_to_base_tag<Derived>::type> && 
      traits<Derived>::direction == 0)
  >
>
  : public Eigen::internal::evaluator_base<Derived>
{
  // Type Defs --------------
  typedef Derived XprType;
  typedef typename Eigen::internal::nested_eval< fornfdm::CSRMatrix, XprType::ColsAtCompileTime >::type ArgTypeNested; // TODO no longer nesting a CSRMatrix...
  typedef typename Eigen::internal::remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 
  enum { CoeffReadCost = Eigen::internal::evaluator<fornfdm::CSRMatrix>::CoeffReadCost, Flags = Eigen::RowMajorBit }; // TODO how important is CoeffReadCost? can I just hard code it? 

  struct InnerIterator
  {
    static constexpr std::size_t N = traits<Derived>::max_arity;
    using Reader = decltype(std::declval<Evaluator<Derived>>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>()));
    // Constructor =============
    InnerIterator(const EigenEvaluatorImpl& eval, Index row_idx)
      : m_eval(eval), 
      m_row_idx(row_idx), 
      m_axis_idx(row_idx % eval.m_xpr.getAxisSize()),
      m_col_offset(eval.m_xpr.getAxisSize() * (row_idx / eval.m_xpr.getAxisSize())), 
      m_counter(0),
      m_offset(eval.m_xpr.getOutersPtr()[m_axis_idx]),
      m_size(eval.m_xpr.getOutersPtr()[m_axis_idx+1] - m_offset),
      m_inner_indices(eval.m_xpr.getInnersPtr() + m_offset),
      m_weights(eval.m_xpr.getWeightsPtr() + m_offset * count_orders<typename traits<Derived>::orders>::value),
      m_reader(eval.m_eval.template createExactReader<N>(fornfdm::Coordinate<N>(eval.m_mesh.get(), row_idx)))
    {}
    // Member Functions ----------------------------------
    operator bool() const { return (m_counter != m_size); }
    void operator++(){ ++m_counter; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_inner_indices[m_counter]; }
    Index index() const { return m_col_offset + m_inner_indices[m_counter]; }
    Scalar value() const { return m_reader(m_weights, m_counter, m_size, typename traits<Derived>::orders{}); }
    // Member Data ----------------------------------------
    const EigenEvaluatorImpl& m_eval;
    std::size_t m_row_idx;
    std::size_t m_axis_idx;
    std::size_t m_counter;
    std::size_t m_col_offset;
    std::size_t m_offset;
    std::size_t m_size;
    const std::size_t* m_inner_indices;
    const fornfdm::Scalar* m_weights;
    Reader m_reader; // depends on coordinate + time
  };

  // Constructors ====================================== 
  EigenEvaluatorImpl(const XprType& xpr)
    : m_xpr(xpr),
    m_eval(xpr, xpr.getTime()),
    m_mesh(xpr.getMesh()),
    m_time(xpr.getTime()),
    m_size(xpr.rows())
  {}

  EigenEvaluatorImpl(const XprType& xpr, fornfdm::Real t)
    : m_xpr(xpr),
    m_eval(xpr, t),
    m_mesh(xpr.getMesh()),
    m_time(t),
    m_size(xpr.rows())
  {}

  // Member Functions -----------------------------------
  Index rows() const {return m_size; }
  Index cols() const {return m_size; }
  Index outerSize() const { return m_size; }
  Index innerSize() const { return m_size; }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Member Data ----------------------------------------
  const Derived& m_xpr;
  fornfdm::Real m_time;
  std::size_t m_size;
  Evaluator<Derived> m_eval;
  SharedConstMesh m_mesh;
};

// With direction != 0 
template<class Derived>
struct EigenEvaluatorImpl<
  Derived, 
  std::enable_if_t<
    (std::is_base_of_v<StoredWeightsTag, typename map_to_base_tag<Derived>::type> && 
      traits<Derived>::direction != 0)
  >
>
  : public Eigen::internal::evaluator_base<Derived>
{
  // Type Defs --------------
  typedef Derived XprType;
  typedef typename Eigen::internal::nested_eval< fornfdm::CSRMatrix, XprType::ColsAtCompileTime >::type ArgTypeNested; // TODO no longer nesting a CSRMatrix...
  typedef typename Eigen::internal::remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 
  enum { CoeffReadCost = Eigen::internal::evaluator<fornfdm::CSRMatrix>::CoeffReadCost, Flags = Eigen::RowMajorBit }; // TODO how important is CoeffReadCost? can I just hard code it? 

  struct InnerIterator
  {
    static constexpr std::size_t N = traits<Derived>::max_arity;
    using Reader = decltype(std::declval<Evaluator<Derived>>().template createExactReader<N>(std::declval<const fornfdm::Coordinate<N>&>()));
    // Constructor =============
    InnerIterator(const EigenEvaluatorImpl& eval, Index row_idx)
      : m_eval(eval), 
      m_row_idx(row_idx), 
      m_axis_idx((row_idx / eval.m_xpr.getProductBefore()) % eval.m_xpr.getAxisSize()),
      m_col_offset(eval.m_xpr.getProductBefore() * eval.m_xpr.getAxisSize() * (row_idx / (eval.m_xpr.getProductBefore() * eval.m_xpr.getAxisSize())) + row_idx % eval.m_xpr.getProductBefore()), 
      m_counter(0),
      m_offset(eval.m_xpr.getOutersPtr()[m_axis_idx]),
      m_size(eval.m_xpr.getOutersPtr()[m_axis_idx+1] - m_offset),
      m_inner_indices(eval.m_xpr.getInnersPtr() + m_offset),
      m_weights(eval.m_xpr.getWeightsPtr() + m_offset * count_orders<typename traits<Derived>::orders>::value),
      m_reader(eval.m_eval.template createExactReader<N>(fornfdm::Coordinate<N>(eval.m_mesh.get(), row_idx)))
    {}
    // Member Functions ----------------------------------
    operator bool() const { return (m_counter != m_size); }
    void operator++(){ ++m_counter; }
    Index row() const { return m_row_idx; }
    Index col() const { return m_col_offset + m_inner_indices[m_counter] * m_eval.m_xpr.getProductBefore(); }
    Index index() const { return m_col_offset + m_inner_indices[m_counter] * m_eval.m_xpr.getProductBefore(); }
    Scalar value() const { return m_reader(m_weights, m_counter, m_size, typename traits<Derived>::orders{}); }
    // Member Data ----------------------------------------
    const EigenEvaluatorImpl& m_eval;
    std::size_t m_row_idx;
    std::size_t m_axis_idx;
    std::size_t m_counter;
    std::size_t m_col_offset;
    std::size_t m_offset;
    std::size_t m_size;
    const std::size_t* m_inner_indices;
    const fornfdm::Scalar* m_weights;
    Reader m_reader; // depends on coordinate + time
  };

  // Constructors ====================================== 
  EigenEvaluatorImpl(const XprType& xpr)
    : m_xpr(xpr),
    m_eval(xpr, xpr.getTime()),
    m_mesh(xpr.getMesh()),
    m_time(xpr.getTime()),
    m_size(xpr.rows())
  {}

  EigenEvaluatorImpl(const XprType& xpr, fornfdm::Real t)
    : m_xpr(xpr),
    m_eval(xpr, t),
    m_mesh(xpr.getMesh()),
    m_time(t),
    m_size(xpr.rows())
  {}

  // Member Functions -----------------------------------
  Index rows() const {return m_size; }
  Index cols() const {return m_size; }
  Index outerSize() const { return m_size; }
  Index innerSize() const { return m_size; }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Member Data ----------------------------------------
  const Derived& m_xpr;
  fornfdm::Real m_time;
  std::size_t m_size;
  Evaluator<Derived> m_eval;
  SharedConstMesh m_mesh;
};

// ==================================================================
// TimeEvaluation<ArgType>
// ArgType --> TimeDepLeftKroneckerTag 
// Use single kronecker product to make block diagonal: I @ D
// ==================================================================

template<class ArgType>
struct EigenEvaluatorImpl<
  fornfdm::linops::TimeEvaluation<ArgType>, 
  std::enable_if_t<
    (std::is_same_v<TimeDepLeftKroneckerTag, typename map_to_base_tag<ArgType>::type> &&
      traits<ArgType>::direction == 0)
  >
> : public Eigen::internal::evaluator_base<fornfdm::linops::TimeEvaluation<ArgType>>
{
  // Type Defs 
  typedef typename Eigen::internal::traits<fornfdm::linops::TimeEvaluation<ArgType>>::StorageIndex StorageIndex; 
  typedef typename Eigen::internal::traits<fornfdm::linops::TimeEvaluation<ArgType>>::Scalar Scalar; 
  typedef typename fornfdm::linops::TimeEvaluation<ArgType> Nested;

  // Flags -------
  enum { CoeffReadCost = Eigen::internal::evaluator<ArgType>::CoeffReadCost, Flags = Eigen::RowMajor };

  // Member Data ---------------- 
  using traits_t = fornfdm::linops::internal::traits<ArgType>;
  fornfdm::linops::internal::Evaluator<ArgType> m_eval;
  Nested m_xpr;
  SharedConstMesh m_mesh;
  StorageIndex m_stencil_size;

  struct InnerIterator
  {
    // Member Data -----------
    const EigenEvaluatorImpl& m_eval; 
    std::size_t m_idx;
    std::size_t m_row_idx;
    std::size_t m_offset; 
    std::size_t m_stride; 
    typename fornfdm::linops::internal::Evaluator<ArgType>::Row m_row;

    // Constructor --------
    InnerIterator(const EigenEvaluatorImpl& eval, std::size_t row_index)
      : m_eval(eval),
      m_idx(0), 
      m_row_idx(row_index),
      m_offset(m_eval.m_xpr.m_arg.getStencil().cols() * (m_row_idx / (m_eval.m_xpr.m_arg.getStencil().rows()))),
      m_row(eval.m_eval, eval.m_mesh.get(), row_index % (m_eval.m_xpr.m_arg.getStencil().rows()), row_index)
    {}
    
    // Member Functions --------------
    operator bool() const { return (m_idx != m_row.size()); }
    void operator++(){ ++m_idx; }
    StorageIndex row() const { return m_row_idx; }
    StorageIndex col() const { return m_offset + m_row.index(m_idx); }
    StorageIndex index() const { return m_offset + m_row.index(m_idx); }
    Scalar value() const { return m_row.value(m_idx); }
  };

  // Constructor ---------------- 
  EigenEvaluatorImpl(const fornfdm::linops::TimeEvaluation<ArgType>& xpr)
    : m_eval(xpr.m_arg, xpr.m_time), m_xpr(xpr), m_mesh(xpr.m_arg.getMesh())
  {}

  // member Functions ------------------ 
  auto rows() const { return m_xpr.rows(); }
  auto cols() const { return m_xpr.cols(); }
  auto outerSize() const { return  m_xpr.rows(); }
  auto innerSize() const { return  m_xpr.cols(); }
  auto nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }
};

// ==================================================================
// TimeEvaluation<ArgType>
// ArgType --> TimeDepDoubleKroneckerTag 
// Use double kronecker product into higher dimension + block diagonal: I @ D @ I
// ==================================================================

template<class ArgType>
struct EigenEvaluatorImpl<
  fornfdm::linops::TimeEvaluation<ArgType>, 
  std::enable_if_t<
    (std::is_same_v<TimeDepDoubleKroneckerTag, typename map_to_base_tag<ArgType>::type> ||
    (std::is_same_v<TimeDepLeftKroneckerTag, typename map_to_base_tag<ArgType>::type> && 
      traits<ArgType>::direction != 0))
  >
> : public Eigen::internal::evaluator_base<fornfdm::linops::TimeEvaluation<ArgType>>
{
  // Type Defs 
  typedef typename Eigen::internal::traits<fornfdm::linops::TimeEvaluation<ArgType>>::StorageIndex StorageIndex; 
  typedef typename Eigen::internal::traits<fornfdm::linops::TimeEvaluation<ArgType>>::Scalar Scalar; 
  typedef typename fornfdm::linops::TimeEvaluation<ArgType> Nested;

  // Flags -------
  enum { CoeffReadCost = Eigen::internal::evaluator<ArgType>::CoeffReadCost, Flags = Eigen::RowMajor };

  // Member Data ---------------- 
  using traits_t = fornfdm::linops::internal::traits<ArgType>;
  fornfdm::linops::internal::Evaluator<ArgType> m_eval;
  Nested m_xpr;
  SharedConstMesh m_mesh;
  StorageIndex m_stencil_size;
  StorageIndex m_prod_before;

  struct InnerIterator
  {
    // Member Data -----------
    const EigenEvaluatorImpl& m_eval; 
    std::size_t m_idx;
    std::size_t m_row_idx;
    std::size_t m_offset; 
    std::size_t m_stride; 
    typename fornfdm::linops::internal::Evaluator<ArgType>::Row m_row;

    // Constructor --------
    InnerIterator(const EigenEvaluatorImpl& eval, std::size_t row_index)
      : m_eval(eval),
      m_idx(0), 
      m_row_idx(row_index),
      m_offset( m_eval.m_prod_before * m_eval.m_stencil_size * (m_row_idx / (m_eval.m_prod_before * m_eval.m_stencil_size)) + (m_row_idx % m_eval.m_prod_before) ),
      m_row(eval.m_eval, eval.m_mesh.get(), (row_index / m_eval.m_prod_before)%(m_eval.m_stencil_size), row_index)
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
  EigenEvaluatorImpl(const fornfdm::linops::TimeEvaluation<ArgType>& xpr)
    : m_eval(xpr.m_arg, xpr.m_time), m_xpr(xpr), m_mesh(xpr.m_arg.getMesh())
  {
    // work around for TimeDepLeftKronecker not exposing a getProductBefore() or correct stencil size
    if constexpr(std::is_same_v<TimeDepDoubleKroneckerTag, typename map_to_base_tag<ArgType>::type>)
    {
      m_stencil_size = xpr.m_arg.getStencil().rows();
      m_prod_before = xpr.m_arg.getProductBefore();
    }
    else
    {
      m_stencil_size = m_mesh->sizeOfDim(traits<ArgType>::direction);
      m_prod_before = m_mesh->sizesMiddleProduct(0, traits<ArgType>::direction);
    }    
  }

  // member Functions ------------------ 
  auto rows() const { return m_xpr.rows(); }
  auto cols() const { return m_xpr.cols(); }
  auto outerSize() const { return  m_xpr.rows(); }
  auto innerSize() const { return  m_xpr.cols(); }
  auto nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }
};

// ==================================================================
// TimeEvaluation<ArgType>
// ArgType --> TimeDepStoredWeightsTag
// ==================================================================

template<class ArgType>
struct EigenEvaluatorImpl<
  fornfdm::linops::TimeEvaluation<ArgType>,
  std::enable_if_t<
    std::is_same_v< TimeDepStoredWeightsTag, typename map_to_base_tag<ArgType>::type >
  >
> : public EigenEvaluatorImpl<ArgType>
{
  // Flags ----- 
  enum { CoeffReadCost = EigenEvaluatorImpl<ArgType>::CoeffReadCost, Flags = Eigen::RowMajor };
  using InnerIterator = typename EigenEvaluatorImpl<ArgType>::InnerIterator;
  EigenEvaluatorImpl(const fornfdm::linops::TimeEvaluation<ArgType>& time_eval)
    : EigenEvaluatorImpl<ArgType>(time_eval.m_arg, time_eval.m_time)
  {}
};

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm 

#endif // EigenEvaluatorImpl.hpp