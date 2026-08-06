// MatrixFree.hpp
//
// decorator class that wraps a linear operator. 
// internall stores the fdm weights after setMesh. 
// doesn't evaluate the into a stencil. 
//
// JAF 8/5/2026

#ifndef FORNFDM_DIFFOPS_MATRIXFREE_H
#define FORNFDM_DIFFOPS_MATRIXFREE_H

#include<type_traits>
#include<memory>
#include "traits.hpp"
#include "EvaluatorBase.hpp"
#include "EigenEvaluatorImpl.hpp"

namespace fornfdm{
  namespace linops{

// forward declare
template<class XprType, class>
class MatrixFree;

namespace internal{

// traits 
template<class XprType>
struct traits_impl<MatrixFree<XprType>>
{
  static constexpr bool is_linop = true; 
  static constexpr bool is_unarop = true; 
  static constexpr bool is_binop = false; 
  static constexpr bool is_ternop = false; 
  static constexpr std::size_t max_arity = traits<XprType>::max_arity; 
  static constexpr bool is_timedep = traits<XprType>::is_timedep; // if either L/R is timedep the xpr is time dep 
  static constexpr int direction = traits<XprType>::direction; // by default mixing operators results in undefined direction... 
  static constexpr std::size_t max_order = traits<XprType>::max_order; // highest order of derivative in the expression 
  typedef typename traits<XprType>::orders orders;
  typedef void node_selector_tag; // give priority to lhs 
};

// Evaluator 
template<class XprType>
struct Evaluator<MatrixFree<XprType>> : public EvaluatorBase<MatrixFree<XprType>>
{
  using XprTypeCleaned = std::remove_cv_t<std::remove_reference_t<XprType>>;
  Evaluator<XprTypeCleaned> m_linop_eval; 

  Evaluator(const MatrixFree<XprType>& xpr, fornfdm::Real t) 
  : EvaluatorBase<MatrixFree<XprType>>(t),
    m_linop_eval(xpr.nestedExpression(), t)
  {}

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord) const
  { 
    return m_linop_eval.template createReader<N>(coord);
  }

  template<std::size_t N>
  auto createExactReader(const fornfdm::Coordinate<N>& coord) const
  {
    return m_linop_eval.template createExactReader<N>(coord);
  }
};

// struct to add a fornfdm::Real conditionally 
template<bool>
struct TimeDepData{};

template<>
struct TimeDepData<true>{ fornfdm::Real m_time; };

} // end namespace internal 
} // end namespace linops
} // end namespace fornfdm

namespace Eigen{
namespace internal{

template<class XprType>
struct traits<fornfdm::linops::MatrixFree<XprType>>
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
    Flags = Eigen::RowMajorBit | NestByRefBit, /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ 
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

template<typename XprType_>
struct evaluator<fornfdm::linops::MatrixFree<XprType_>> 
  : public fornfdm::linops::internal::EigenEvaluatorImpl<fornfdm::linops::MatrixFree<XprType_>>
{
  using XprType = fornfdm::linops::MatrixFree<XprType_>; 
  using Impl = typename fornfdm::linops::internal::EigenEvaluatorImpl<XprType>; 
  enum {CoeffReadCost = Impl::CoeffReadCost, Flags = Impl::Flags};
  using InnerIterator = typename Impl::InnerIterator; 
  evaluator(const XprType& xpr)
    : Impl(xpr)
  {}
}; 

} // end namespac internal 
} // end namespac Eigen 

namespace fornfdm{
namespace linops{

// ==================================================================
// MatrixFree (direction == 0)
// ==================================================================

template<class XprType>
class MatrixFree<
  XprType, 
  std::enable_if_t<
    (internal::traits<XprType>::direction == 0)
  >
> : public Eigen::SparseMatrixBase<MatrixFree<XprType>>, 
private internal::TimeDepData<internal::traits<XprType>::is_timedep>
{
  public:
    // Type Defs -----------------
    typedef Eigen::SparseMatrixBase<MatrixFree<XprType>> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(MatrixFree) 
    using XprTypeNested = typename fornfdm::linops::internal::NestedStorage<XprType>::type;
    using NestedExpression = typename std::remove_reference<std::remove_cv_t<XprType>>::type;

  protected:
    // Member Data ---------------
    XprTypeNested m_nested;
    WeakConstMesh m_mesh_observed;
    std::size_t m_axis_size=0;
    std::size_t m_nnz=0;
    std::size_t m_prod_after;
    std::unique_ptr<std::size_t[]> m_outer_ptr=nullptr;
    std::unique_ptr<std::size_t[]> m_inner_ptr=nullptr;
    std::unique_ptr<fornfdm::Scalar[]> m_weights_ptr=nullptr;

  public:
    // Constructors + Destructor 
    MatrixFree()=delete;
    template<class Xpr, typename = std::enable_if_t<internal::is_partialderiv_crtp<Xpr>::value>>
    MatrixFree(Xpr&& linop)
      : m_nested(std::forward<Xpr>(linop))
    {}
    MatrixFree(const MatrixFree& other)
      : m_nested(other.m_nested),
    m_mesh_observed(other.m_mesh_observed),
    m_axis_size(other.m_axis_size),
    m_nnz(other.m_nnz),
    m_prod_after(other.m_prod_after)
    {
      // unique_ptrs must deep copy.
      if(other.m_outer_ptr){
        m_outer_ptr.reset(new std::size_t[ m_axis_size + 1 ]);
        for(auto i=0; i<(m_axis_size+1); ++i){
          m_outer_ptr.get()[i] = other.m_outer_ptr.get()[i];
        }
      }
      if(other.m_inner_ptr){
        m_inner_ptr.reset(new std::size_t[ m_nnz ]);
        for(auto i=0; i<(m_nnz); ++i){
          m_inner_ptr.get()[i] = other.m_inner_ptr.get()[i];
        }
      }
      if(other.m_weights_ptr){
        constexpr std::size_t count_of_orders = internal::count_orders< typename internal::traits<XprType>::orders>::value;
        m_weights_ptr.reset(new fornfdm::Scalar[ m_nnz * count_of_orders]);
        for(auto i=0; i<(m_nnz * count_of_orders); ++i){
          m_weights_ptr.get()[i] = other.m_weights_ptr.get()[i];
        }
      }
    } 
    MatrixFree(MatrixFree&& other)
      : m_nested(other.m_nested),
      m_mesh_observed(other.m_mesh_observed),
      m_axis_size(other.m_axis_size),
      m_nnz(other.m_nnz),
      m_prod_after(other.m_prod_after),
      m_outer_ptr(std::move(other.m_outer_ptr)),  
      m_inner_ptr(std::move(other.m_inner_ptr)),  
      m_weights_ptr(std::move(other.m_weights_ptr))  
    {}
    ~MatrixFree()=default;

    // Member Funcs ---------
    const auto& nestedExpression() const { return m_nested; }
    auto rows() const { return m_axis_size * m_prod_after; }
    auto cols() const { return m_axis_size * m_prod_after; }
    auto nonZerosEstimate() const { return m_prod_after * m_nnz; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<MatrixFree<XprType>>*>(this); }

    auto getProductAfter() const { return m_prod_after; }
    auto getAxisSize() const { return m_axis_size; }

    inline const auto* getOutersPtr() const { return m_outer_ptr.get(); } // TODO best to not store in memory ??? NodeSelector maps size + index to inner offset
    inline const auto* getInnersPtr() const { return m_inner_ptr.get(); }
    inline const auto* getWeightsPtr() const { return m_weights_ptr.get(); }

    fornfdm::Real getTime() const 
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        return this->m_time;
      } 
      else{
        return -1.0;
      }
    }

    const auto& evalTime(fornfdm::Real t) const & 
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        return fornfdm::linops::TimeEvaluation(*this, t);
      } 
      else{
        // nothing to calculate!
        return toEigen(); 
      }
    } 

    void setTime(fornfdm::Real t)
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        this->m_time = t;
      } 
    }

    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }

    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<XprType>;
      using Selector = internal::NodeSelector< typename traits_t::node_selector_tag, traits_t::max_order+1>;
      constexpr std::size_t count_of_orders = internal::count_orders< typename traits_t::orders >::value;
      m_mesh_observed = m;
      m_prod_after = m->sizesMiddleProduct(traits_t::direction + 1, m->numDims());
      const fornfdm::Vector& axis = m->getAxis(traits_t::direction);
      std::size_t axis_size = m->sizeOfDim(traits_t::direction);
      std::size_t nnz = Selector::sumNodesPerRow(axis);

      if(axis_size > m_axis_size){
        m_outer_ptr.reset(new std::size_t[ axis_size + 1 ]);
      }
      m_axis_size = axis_size;

      if(nnz > m_nnz){
        m_inner_ptr.reset(new std::size_t[ nnz ]);
        m_weights_ptr.reset(new fornfdm::Scalar[ nnz * count_of_orders]);
      }
      m_nnz = nnz; 

      for(auto row_idx = 0; row_idx < axis_size; ++row_idx)
      {
        Selector nodes(axis, row_idx);
        
        // write the outer  ptr 
        m_outer_ptr[row_idx] = nodes.nonZerosOffset;
        
        // write the inners 
        std::copy(nodes.nodeIndices.cbegin(), nodes.nodeIndices.cend(), std::next(m_inner_ptr.get(),nodes.nonZerosOffset));
        
        // write the weights.
        using Calc = fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::max_order>; 
        Calc calc(nodes.x_bar, nodes.nodeValues.cbegin(), nodes.nodeValues.cend());
        assignment_helper(
          m_weights_ptr.get() + nodes.nonZerosOffset * count_of_orders, 
          calc.getArray().data(),
          calc.getNumNodesUsed(), 
          typename traits_t::orders{}
        );
      }
      m_outer_ptr[axis_size] = nnz;
    }

  private:
    template<std::size_t... Idxs>
    void assignment_helper(fornfdm::Scalar* dest, const fornfdm::Scalar* src, std::size_t n, std::index_sequence<Idxs...>)
    {
      ((dest = std::copy(src + (n*Idxs), src + (n*(Idxs+1)), dest)), ...);
    }
};

// ==================================================================
// MatrixFree (direction != 0)
// ==================================================================

template<class XprType>
class MatrixFree<
  XprType, 
  std::enable_if_t<
    (internal::traits<XprType>::direction != 0)
  >
> : public Eigen::SparseMatrixBase<MatrixFree<XprType>>, 
private internal::TimeDepData<internal::traits<XprType>::is_timedep>
{
  public:
    // Type Defs -----------------
    typedef Eigen::SparseMatrixBase<MatrixFree<XprType>> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(MatrixFree) 
    using XprTypeNested = typename fornfdm::linops::internal::NestedStorage<XprType>::type;
    using NestedExpression = typename std::remove_reference<std::remove_cv_t<XprType>>::type;

  protected:
    // Member Data ---------------
    XprTypeNested m_nested;
    WeakConstMesh m_mesh_observed;
    std::size_t m_axis_size=0;
    std::size_t m_nnz=0;
    std::size_t m_prod_before;
    std::size_t m_prod_after;
    std::unique_ptr<std::size_t[]> m_outer_ptr=nullptr;
    std::unique_ptr<std::size_t[]> m_inner_ptr=nullptr;
    std::unique_ptr<fornfdm::Scalar[]> m_weights_ptr=nullptr;

  public:
    // Constructors + Destructor 
    MatrixFree()=delete;
    template<class Xpr, typename = std::enable_if_t<internal::is_partialderiv_crtp<Xpr>::value>>
    MatrixFree(Xpr&& linop)
      : m_nested(std::forward<Xpr>(linop))
    {}
    MatrixFree(const MatrixFree& other)
      : m_nested(other.m_nested),
    m_mesh_observed(other.m_mesh_observed),
    m_axis_size(other.m_axis_size),
    m_nnz(other.m_nnz),
    m_prod_before(other.m_prod_before),
    m_prod_after(other.m_prod_after)
    {
      // unique_ptrs must deep copy.
      if(other.m_outer_ptr){
        m_outer_ptr.reset(new std::size_t[ m_axis_size + 1 ]);
        for(auto i=0; i<(m_axis_size+1); ++i){
          m_outer_ptr.get()[i] = other.m_outer_ptr.get()[i];
        }
      }
      if(other.m_inner_ptr){
        m_inner_ptr.reset(new std::size_t[ m_nnz ]);
        for(auto i=0; i<(m_nnz); ++i){
          m_inner_ptr.get()[i] = other.m_inner_ptr.get()[i];
        }
      }
      if(other.m_weights_ptr){
        constexpr std::size_t count_of_orders = internal::count_orders< typename internal::traits<XprType>::orders>::value;
        m_weights_ptr.reset(new fornfdm::Scalar[ m_nnz * count_of_orders]);
        for(auto i=0; i<(m_nnz * count_of_orders); ++i){
          m_weights_ptr.get()[i] = other.m_weights_ptr.get()[i];
        }
      }
    } 
    MatrixFree(MatrixFree&& other)
      : m_nested(other.m_nested),
      m_mesh_observed(other.m_mesh_observed),
      m_axis_size(other.m_axis_size),
      m_nnz(other.m_nnz),
      m_prod_before(other.m_prod_before),
      m_prod_after(other.m_prod_after),
      m_outer_ptr(std::move(other.m_outer_ptr)),  
      m_inner_ptr(std::move(other.m_inner_ptr)),  
      m_weights_ptr(std::move(other.m_weights_ptr))
    {}
    ~MatrixFree()=default;

    // Member Funcs ---------
    const auto& nestedExpression() const { return m_nested; }
    auto rows() const { return m_prod_before * m_axis_size * m_prod_after; }
    auto cols() const { return m_prod_before * m_axis_size * m_prod_after; }
    auto nonZerosEstimate() const { return m_prod_after * m_nnz * m_prod_before; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<MatrixFree<XprType>>*>(this); }

    auto getProductBefore() const { return m_prod_before; }
    auto getProductAfter() const { return m_prod_after; }
    auto getAxisSize() const { return m_axis_size; }

    inline const auto* getOutersPtr() const { return m_outer_ptr.get(); } // TODO best to not store in memory ??? NodeSelector maps size + index to inner offset
    inline const auto* getInnersPtr() const { return m_inner_ptr.get(); }
    inline const auto* getWeightsPtr() const { return m_weights_ptr.get(); }

    fornfdm::Real getTime() const 
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        return this->m_time;
      } 
      else{
        return -1.0;
      }
    }

    const auto& evalTime(fornfdm::Real t) const & 
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        return fornfdm::linops::TimeEvaluation(*this, t);
      } 
      else{
        // nothing to calculate!
        return toEigen(); 
      }
    } 

    void setTime(fornfdm::Real t)
    {
      if constexpr(internal::traits<XprType>::is_timedep){
        this->m_time = t;
      } 
    }

    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }

    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<XprType>;
      using Selector = internal::NodeSelector< typename traits_t::node_selector_tag, traits_t::max_order+1>;
      constexpr std::size_t count_of_orders = internal::count_orders< typename traits_t::orders >::value;
      m_mesh_observed = m;
      m_prod_before = m->sizesMiddleProduct(0,traits_t::direction); 
      m_prod_after = m->sizesMiddleProduct(traits_t::direction + 1, m->numDims());
      const fornfdm::Vector& axis = m->getAxis(traits_t::direction);
      std::size_t axis_size = m->sizeOfDim(traits_t::direction);
      std::size_t nnz = Selector::sumNodesPerRow(axis);

      if(axis_size > m_axis_size){
        m_outer_ptr.reset(new std::size_t[ axis_size + 1 ]);
      }
      m_axis_size = axis_size;

      if(nnz > m_nnz){
        m_inner_ptr.reset(new std::size_t[ nnz ]);
        m_weights_ptr.reset(new fornfdm::Scalar[ nnz * count_of_orders]);
      }
      m_nnz = nnz; 

      for(auto row_idx = 0; row_idx < axis_size; ++row_idx)
      {
        Selector nodes(axis, row_idx);
        
        // write the outer  ptr 
        m_outer_ptr[row_idx] = nodes.nonZerosOffset;
        
        // write the inners 
        std::copy(nodes.nodeIndices.cbegin(), nodes.nodeIndices.cend(), std::next(m_inner_ptr.get(),nodes.nonZerosOffset));
        
        // write the weights.
        using Calc = fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::max_order>; 
        Calc calc(nodes.x_bar, nodes.nodeValues.cbegin(), nodes.nodeValues.cend());
        assignment_helper(
          m_weights_ptr.get() + nodes.nonZerosOffset * count_of_orders, 
          calc.getArray().data(),
          calc.getNumNodesUsed(), 
          typename traits_t::orders{}
        );
      }
      m_outer_ptr[axis_size] = nnz;
    }

  private:
    template<std::size_t... Idxs>
    void assignment_helper(fornfdm::Scalar* dest, const fornfdm::Scalar* src, std::size_t n, std::index_sequence<Idxs...>)
    {
      ((dest = std::copy(src + (n*Idxs), src + (n*(Idxs+1)), dest)), ...);
    }
};

// CTAD Guideline
template<class Xpr>
MatrixFree(Xpr&&) -> MatrixFree<Xpr>;

} // end namespace linops 
} // end namespace fornfdm

#endif // MatrixFree.hpp