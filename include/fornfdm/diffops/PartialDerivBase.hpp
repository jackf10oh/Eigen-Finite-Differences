// PartialDerivBase.hpp
//
// CRTP base class all Partial Derivatives and expressions will derive from.
// Expressions are only able to be constructed if they add, minus, scalar mult, etc...
// on the same "selected nodes" in space. (Nwise...Op ~ Node Wise binar/unar Operation)
// and are the same direction partial derivative.  
//
// JAF 4/17/2026 

#ifndef FORNFDM_DIFFOPS_PARTIALDERIVBASE_H
#define FORNFDM_DIFFOPS_PARTIALDERIVBASE_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include "../types.hpp"
#include "traits.hpp"
#include "EvaluatorBase.hpp"

namespace fornfdm{
namespace linops{
namespace internal{

// Evaluator 
template<class Derived>
struct Evaluator<PartialDerivBase<Derived>> : public EvaluatorBase<PartialDerivBase<Derived>>
{
  Evaluator<Derived> m_derived_eval; 

  Evaluator(const PartialDerivBase<Derived>& xpr) : m_derived_eval(xpr.derived()){}

  template<std::size_t N>
  auto createReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  { 
    return m_derived_eval.template createReader<N>(coord,t);
  }

  template<std::size_t N>
  auto createExactReader(const fornfdm::Coordinate<N>& coord, fornfdm::Real t) const
  {
    return m_derived_eval.template createExactReader<N>(coord,t);
  }
}; 

// linops traits
template<class Derived>
struct traits_impl<fornfdm::linops::PartialDerivBase<Derived>> : traits_impl<Derived>{}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm

namespace Eigen{
namespace internal{

// Eigen traits of PartialDerivBase is same as Derived
template<class Derived>
struct traits<fornfdm::linops::PartialDerivBase<Derived>> : public traits<Derived>{}; 

// Eigen evaluator of PartialDerivBase is same as Derived.
template<class Derived>
struct evaluator<fornfdm::linops::PartialDerivBase<Derived>> 
  : public evaluator<Derived>
{
  using InnerIterator = typename evaluator<Derived>::InnerIterator; 
  evaluator(const fornfdm::linops::PartialDerivBase<Derived>& xpr)
    : evaluator<Derived>(xpr.derived())
  {}
}; 

} // end namespace internal 
} // end namespace Eigen

namespace fornfdm{
namespace linops{

// empty base class. 
template<class Derived, class TagType>
class PartialDerivBase{}; 

// ==================================================================
// map_to_base_tag isn't specialized. everything delegates to derived!
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived,void> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

    // Member Functions ----------------------------- 
    void setMesh(const std::shared_ptr<const Mesh>& m){ derived().setMesh(m); }
    SharedConstMesh getMesh() const { return derived().getMesh(); }
    decltype(auto) evalTime(fornfdm::Real t) const & { return derived().evalTime(t); }
    void setTime(fornfdm::Real t){ derived().setTime(t); }
    fornfdm::Real getTime() const { return derived().getTime(); }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }
};

// ==================================================================
// LeftKroneckerTag
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::LeftKroneckerTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:
    // Member Data ---------------- 
    std::size_t m_prod_after; 
    WeakConstMesh m_mesh_observed;
    fornfdm::CSRMatrix m_stencil;

  public:
    // Member Functions 
    void setMesh(const SharedConstMesh& m)
    {
      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      assert((0 <  m->numDims()) && "direction must be < # of dims.");
      
      // store a weak_ptr to the mesh
      m_mesh_observed = m;

      // produce m_product_after repeats of m_stencil on the diagonal
      m_prod_after = m->sizesMiddleProduct(1,m->numDims()); 

      // setup hot loop
      using Evaluator = fornfdm::linops::internal::Evaluator<Derived>; 
      Evaluator eval(derived()); 
      const auto& axis = m->getAxis(0); 
      const std::size_t axis_size = m->sizeOfDim(0);  

      // resize + reserve the stencil 
      std::size_t nnz = Evaluator::totalNonZeros(axis); 
      m_stencil.reserve(nnz); 
      m_stencil.resize(axis_size, axis_size);
        
      // write node wise expressions into each row stencil 
      for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
      {
        // use node selector 
        typename Evaluator::Row row(eval, m.get(), row_idx, row_idx, -1.0); 
        // copy the indices into m_stencil's inner indices ptr
        m_stencil.outerIndexPtr()[row_idx] = row.offset(); 
        for(auto i=0; i<row.size(); ++i){
          m_stencil.innerIndexPtr()[row.offset() + i] = row.index(i);
          m_stencil.valuePtr()[row.offset() + i] = row.value(i);
        }
      }
      // very last outer index needs set. 
      m_stencil.outerIndexPtr()[axis_size] = nnz; 
    }
    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    const auto& evalTime(fornfdm::Real) const & { return toEigen(); } // nothing to calculate!
    void setTime(fornfdm::Real){/*non op*/}
    fornfdm::Real getTime() const { return -1.0; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }
    auto rows() const { return m_prod_after * m_stencil.rows(); }
    auto cols() const { return m_prod_after * m_stencil.cols(); }
    auto nonZerosEstimate() const { return m_prod_after * m_stencil.nonZeros(); }
    const fornfdm::CSRMatrix& getStencil() const { return m_stencil; }
    std::size_t getProductAfter() const { return m_prod_after; }
};

// ==================================================================
// TimeDepLeftKroneckerTag
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::TimeDepLeftKroneckerTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:
    // Member Data ---------------- 
    std::size_t m_prod_after; 
    WeakConstMesh m_mesh_observed;
    fornfdm::CSRMatrix m_stencil;
    const fornfdm::Mesh* m_mesh_raw;
    fornfdm::Real m_time = -1.0; 
    bool m_need_col_indices; // signals the first time setTime() is called after a new mesh

  public:
    // Member Functions 
    void setMesh(const SharedConstMesh& m)
    {
      m_mesh_observed = m;
      m_mesh_raw = m.get();
      // produce m_product_after repeats of m_stencil on the diagonal
      m_prod_after = m->sizesMiddleProduct(1,m->numDims());
      // setup
      const std::size_t axis_size = m->sizeOfDim(0);
      // resize + reserve the stencil
      m_stencil.resize(axis_size, axis_size);
      m_need_col_indices = true;
    }
    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    auto evalTime(fornfdm::Real t) const & { return fornfdm::linops::TimeEvaluation<Derived>(derived(), t); }
    void setTime(fornfdm::Real t)
    {
      // store the new time
      m_time = t;

      // setup hot loop
      using Evaluator = fornfdm::linops::internal::Evaluator<Derived>; 
      Evaluator eval(derived()); 
      const fornfdm::Vector& axis = m_mesh_raw->getAxis(0);
      std::size_t axis_size = m_stencil.rows(); 

      if(m_need_col_indices)
      {
        // assignment loop. with assignment to innerIndexPtr
        std::size_t nnz = Evaluator::totalNonZeros(axis); 
        m_stencil.reserve(nnz);

        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          // use node selector 
          typename Evaluator::Row row(eval, m_mesh_raw, row_idx, row_idx, t);
          // copy the indices into m_stencil's inner indices ptr
          m_stencil.outerIndexPtr()[row_idx] = row.offset(); 
          for(auto i=0; i<row.size(); ++i){
            m_stencil.innerIndexPtr()[row.offset() + i] = row.index(i);
            m_stencil.valuePtr()[row.offset() + i] = row.value(i);
          }
        }
        // very last outer index needs set. 
        m_stencil.outerIndexPtr()[axis_size] = nnz; 
        m_need_col_indices=false;
      }
      else
      {
        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          typename Evaluator::Row row(eval, m_mesh_raw, row_idx, row_idx, t);
          for(auto i=0; i<row.size(); ++i){
            m_stencil.valuePtr()[row.offset() + i] = row.value(i);
          }
        }
      }
    }
    fornfdm::Real getTime() const { return m_time; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }
    auto rows() const { return m_prod_after * m_stencil.rows(); }
    auto cols() const { return m_prod_after * m_stencil.cols(); }
    const fornfdm::CSRMatrix& getStencil() const { return m_stencil; }
    auto nonZerosEstimate() const { return m_prod_after * m_stencil.nonZeros(); }
    std::size_t getProductAfter() const { return m_prod_after; }
};

// ==================================================================
// DoubleKroneckerTag 
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::DoubleKroneckerTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:
    // Member Data ---------------- 
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 
    WeakConstMesh m_mesh_observed;
    fornfdm::CSRMatrix m_stencil;

  public:
    // Member Functions 
    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<Derived>; 
      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      assert((traits_t::direction <  m->numDims()) && "direction must be < # of dims.");
      
      // store a weak_ptr to the mesh
      m_mesh_observed = m;

      // produce m_product_after repeats of m_stencil on the diagonal
      m_prod_before = m->sizesMiddleProduct(0,traits_t::direction); 
      m_prod_after = m->sizesMiddleProduct(traits_t::direction + 1, m->numDims());

      // setup hot loop
      using Evaluator = fornfdm::linops::internal::Evaluator<Derived>; 
      Evaluator eval(derived()); 
      const auto& axis = m->getAxis(traits_t::direction); 
      const std::size_t axis_size = m->sizeOfDim(traits_t::direction);  

      // resize + reserve the stencil 
      std::size_t nnz = Evaluator::totalNonZeros(axis); 
      m_stencil.reserve(nnz); 
      m_stencil.resize(axis_size, axis_size);
        
      // write node wise expressions into each row stencil 
      for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
      {
        // use node selector 
        typename Evaluator::Row row(eval, m.get(), row_idx, row_idx, -1.0); 
        // copy the indices into m_stencil's inner indices ptr
        m_stencil.outerIndexPtr()[row_idx] = row.offset(); 
        for(auto i=0; i<row.size(); ++i){
          m_stencil.innerIndexPtr()[row.offset() + i] = row.index(i);
          m_stencil.valuePtr()[row.offset() + i] = row.value(i);
        }
      }
      // very last outer index needs set. 
      m_stencil.outerIndexPtr()[axis_size] = nnz; 
    }
    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    const auto& evalTime() const & { return toEigen(); } // nothing to calculate!
    void setTime(fornfdm::Real){/*non op*/}
    fornfdm::Real getTime() const { return -1.0; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }
    auto rows() const { return m_prod_after * m_stencil.rows() * m_prod_before; }
    auto cols() const { return m_prod_after * m_stencil.cols() * m_prod_before; }
    auto nonZerosEstimate() const { return m_prod_after * m_stencil.nonZeros() * m_prod_before; }
    const fornfdm::CSRMatrix& getStencil() const { return m_stencil; }
    std::size_t getProductBefore() const { return m_prod_before; }
    std::size_t getProductAfter() const { return m_prod_after; }
};

// ==================================================================
// TimeDepDoubleKroneckerTag 
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::TimeDepDoubleKroneckerTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:
    // Member Data ---------------- 
    std::size_t m_prod_before; 
    std::size_t m_prod_after;
    WeakConstMesh m_mesh_observed;
    fornfdm::CSRMatrix m_stencil;
    const fornfdm::Mesh* m_mesh_raw;
    fornfdm::Real m_time = -1.0;
    bool m_need_col_indices; // signals the first time setTime() is called after a new mesh

  public:
    // Member Functions 
    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<Derived>; 
      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      assert((traits_t::direction <  m->numDims()) && "direction must be < # of dims.");
      
      m_mesh_observed = m;
      m_mesh_raw = m.get();
      // produce m_product_after repeats of m_stencil on the diagonal
      m_prod_before = m->sizesMiddleProduct(0,traits_t::direction); 
      m_prod_after = m->sizesMiddleProduct(traits_t::direction + 1, m->numDims());

      // resize the stencil
      const std::size_t axis_size = m->sizeOfDim(traits_t::direction);
      m_stencil.resize(axis_size, axis_size);
      m_need_col_indices = true;
    }
    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    auto evalTime(fornfdm::Real t) const & { return fornfdm::linops::TimeEvaluation<Derived>(derived(), t); }
    void setTime(fornfdm::Real t)
    {
      // store the new time
      m_time = t;

      // setup hot loop
      using traits_t = fornfdm::linops::internal::traits<Derived>; 
      using Evaluator = fornfdm::linops::internal::Evaluator<Derived>; 
      Evaluator eval(derived()); 
      const fornfdm::Vector& axis = m_mesh_raw->getAxis(traits_t::direction);
      std::size_t axis_size = m_stencil.rows(); 

      if(m_need_col_indices)
      {
        // assignment loop. with assignment to innerIndexPtr
        std::size_t nnz = Evaluator::totalNonZeros(axis); 
        m_stencil.reserve(nnz);

        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          // use node selector 
          typename Evaluator::Row row(eval, m_mesh_raw, row_idx, row_idx, t);
          // copy the indices into m_stencil's inner indices ptr
          m_stencil.outerIndexPtr()[row_idx] = row.offset(); 
          for(auto i=0; i<row.size(); ++i){
            m_stencil.innerIndexPtr()[row.offset() + i] = row.index(i);
            m_stencil.valuePtr()[row.offset() + i] = row.value(i);
          }
        }
        // very last outer index needs set. 
        m_stencil.outerIndexPtr()[axis_size] = nnz; 
        m_need_col_indices=false;
      }
      else
      {
        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          typename Evaluator::Row row(eval, m_mesh_raw, row_idx, row_idx, t);
          for(auto i=0; i<row.size(); ++i){
            m_stencil.valuePtr()[row.offset() + i] = row.value(i);
          }
        }
      }
    }
    fornfdm::Real getTime() const { return m_time; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }
    auto rows() const { return m_prod_after * m_stencil.rows() * m_prod_before; }
    auto cols() const { return m_prod_after * m_stencil.cols() * m_prod_before; }
    auto nonZerosEstimate() const { return m_prod_after * m_stencil.nonZeros() * m_prod_before; }
    const fornfdm::CSRMatrix& getStencil() const { return m_stencil; }
    std::size_t getProductBefore() const { return m_prod_before; }
    std::size_t getProductAfter() const { return m_prod_after; }
};

// ==================================================================
// StoredWeightsTag 
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::StoredWeightsTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs -----------------
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:

    // Member Data ---------------
    WeakConstMesh m_mesh_observed;
    std::size_t m_axis_size=0;
    std::size_t m_nnz=0;
    std::size_t m_prod_before;
    std::size_t m_prod_after;
    std::unique_ptr<std::size_t[]> m_outer_ptr=nullptr;
    std::unique_ptr<std::size_t[]> m_inner_ptr=nullptr;
    std::unique_ptr<fornfdm::Scalar[]> m_weights_ptr=nullptr;

  public:
    // Member Funcs ---------
    auto rows() const { return m_prod_before * m_axis_size * m_prod_after;};
    auto cols() const { return m_prod_before * m_axis_size * m_prod_after;};
    auto nonZerosEstimate() const { return m_prod_after * m_nnz * m_prod_before; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }

    auto getProductBefore() const { return m_prod_before; }
    auto getProductAfter() const { return m_prod_after; }
    auto getAxisSize() const { return m_axis_size; }

    inline const auto* getOutersPtr() const { return m_outer_ptr.get(); } // TODO best to not store in memory ??? NodeSelector maps size + index to inner offset
    inline const auto* getInnersPtr() const { return m_inner_ptr.get(); }
    inline const auto* getWeightsPtr() const { return m_weights_ptr.get(); }

    fornfdm::Real getTime() const { return -1.0; }
    const auto& evalTime(fornfdm::Real) const & {return toEigen(); } // nothing to calculate!
    void setTime(fornfdm::Real){/*non op*/}

    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<Derived>;
      using Selector = internal::NodeSelector< typename traits_t::node_selector_tag, traits_t::maxOrder+1>;
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
        using Calc = fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::maxOrder>; 
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
      (std::copy(src + (n*Idxs), src + (n*(Idxs+1)), dest), ...);
    }
};

// ==================================================================
// TimeDepStoredWeightsTag
// ==================================================================

template<class Derived>
class PartialDerivBase<Derived, internal::TimeDepStoredWeightsTag> : public Eigen::SparseMatrixBase<Derived>
{
  public:
    // Type Defs -----------------
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

  protected:

    // Member Data ---------------
    WeakConstMesh m_mesh_observed;
    std::size_t m_axis_size=0;
    std::size_t m_nnz=0;
    std::size_t m_prod_before;
    std::size_t m_prod_after;
    std::unique_ptr<std::size_t[]> m_outer_ptr=nullptr;
    std::unique_ptr<std::size_t[]> m_inner_ptr=nullptr;
    std::unique_ptr<fornfdm::Scalar[]> m_weights_ptr=nullptr;
    fornfdm::Real m_time = -1.0;

  public:
    // Member Funcs ---------
    auto rows() const { return m_prod_before * m_axis_size * m_prod_after;};
    auto cols() const { return m_prod_before * m_axis_size * m_prod_after;};
    auto nonZerosEstimate() const { return m_prod_after * m_nnz * m_prod_before; }
    const auto& toEigen() const { return *static_cast<const Eigen::SparseMatrixBase<Derived>*>(this); }

    auto getProductBefore() const { return m_prod_before; }
    auto getProductAfter() const { return m_prod_after; }
    auto getAxisSize() const { return m_axis_size; }

    inline const auto* getOutersPtr() const { return m_outer_ptr.get(); } // TODO best to not store in memory ??? NodeSelector maps size + index to inner offset
    inline const auto* getInnersPtr() const { return m_inner_ptr.get(); }
    inline const auto* getWeightsPtr() const { return m_weights_ptr.get(); }

    fornfdm::Real getTime() const { return m_time; }
    auto evalTime(fornfdm::Real t) const & { return fornfdm::linops::TimeEvaluation<Derived>(derived(), t); }
    void setTime(fornfdm::Real t){ m_time = t; }

    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    void setMesh(const SharedConstMesh& m)
    {
      using traits_t = fornfdm::linops::internal::traits<Derived>;
      using Selector = internal::NodeSelector< typename traits_t::node_selector_tag, traits_t::maxOrder+1>;
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
        using Calc = fornfdm::utils::FornbergStackCalc<Selector::numNodesMax, traits_t::maxOrder>; 
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
      (std::copy(src + (n*Idxs), src + (n*(Idxs+1)), dest), ...);
    }
};

} // end namespace linops
} // end namespace fornfdm 

#endif // PartialDerivBase.hpp  