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

#include<Eigen/SparseCore>
#include "EvaluatorBase.hpp"
#include "../Types.hpp"

namespace fornfdm{
namespace linops{

template<class Derived> class PartialDerivBase; 

namespace internal{ 

// Member data that is only used by time dependent partial derivs 
template<bool isTimeDep>
struct TimeDepMemberData{}; 

template<>
struct TimeDepMemberData<true>{
  fornfdm::Real m_current_time; 
  const Mesh* m_mesh_raw; 

}; 

// Evaluator 
template<class Derived>
struct Evaluator<PartialDerivBase<Derived>> : public EvaluatorBase<PartialDerivBase<Derived>>
{
  Evaluator<Derived> m_derived_eval; 

  Evaluator(const PartialDerivBase<Derived>& xpr) : m_derived_eval(xpr.derived()){}

  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fornfdm::Scalar* weights, std::size_t weights_per_order, const fornfdm::Coordinate<N>& coords)
  {
    return m_derived_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords); 
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

} // end namespace internal 
} // end namespace Eigen

namespace fornfdm{
namespace linops{

template<class Derived>
class PartialDerivBase : public Eigen::SparseMatrixBase<Derived>, protected linops::internal::TimeDepMemberData<linops::internal::traits<Derived>::is_timedep>
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

    // Friends ------------------- 
    friend Eigen::internal::evaluator<PartialDerivBase>; 
    friend fornfdm::linops::internal::EvaluatorBase<PartialDerivBase>; 
    friend fornfdm::linops::internal::Evaluator<PartialDerivBase>; 

  public: // TODO make protected 
    // Member Data ----------------------------------------------
    // const Mesh* m_mesh_raw = nullptr; // This is unused until I need time dependent operators....... 
    std::weak_ptr<const Mesh> m_mesh_observed = {/*nullptr*/}; 
    fornfdm::CSRMatrix m_stencil = {}; 
    std::size_t m_prod_before = 1; 
    std::size_t m_prod_after = 1; 

  public:
    // Member Functions =====================================================
    // FDM Interface -------
    void setMesh(const std::shared_ptr<const Mesh>& m)
    {
      if constexpr(fornfdm::linops::internal::traits<Derived>::is_timedep){
        // just change the mesh thats observed
        this->m_mesh_observed = m; 
        this->m_mesh_raw = m.get(); 
      }
      else{
        // change mesh observed + update m_stencil matrix 
        this->m_mesh_observed = m; 
        setMesh_impl(m.get()); 
      }
    }
    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    void setTime(fornfdm::Real t){ 
      if constexpr(fornfdm::linops::internal::traits<Derived>::is_timedep){
        // update m_current_time + update m_stencil matrix 
        this->m_current_time=t; 
        derived().setTime_hooked(t); 
        setMesh_impl(this->m_mesh_raw); 
      }
    }
    void setTime_hooked(fornfdm::Real t){}
    fornfdm::Real getTime() const {
      if constexpr(fornfdm::linops::internal::traits<Derived>::is_timedep){
        return this->m_current_time; 
      }
      return -1.0; 
    }
    
    // Eigen Interface ------- 
    const auto& toEigen() const { return *static_cast<const Base*>(this); } // prevents custom operators from fornfdm library taking effect. 
    StorageIndex rows() const { return m_prod_before * m_prod_after * m_stencil.rows(); }
    StorageIndex cols() const { return m_prod_before * m_prod_after * m_stencil.cols(); }
    StorageIndex nonZerosEstimate() const {return m_prod_before * m_prod_after * m_stencil.nonZeros(); }

    // Operators ====================================================================== 
    
    // removed assignment from inheritance hierarchy 
    template<typename OtherDerived>
    Derived& operator=(const Eigen::EigenBase<OtherDerived> &other)=delete;

    template<typename OtherDerived>
    inline Derived& assign(const OtherDerived& other)=delete;

    template<typename OtherDerived>
    inline Derived& assignGeneric(const OtherDerived& other)=delete;

  protected:
    // Implementations ----------------------------------------------------------------
    void setMesh_impl(const Mesh* m)
    {
      using traits_t = fornfdm::linops::internal::traits<Derived>; 
      const auto& axis = m->getAxis(traits_t::direction); 
      const std::size_t axis_size = m->sizeOfDim(traits_t::direction); 

      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      bool callable_check = traits_t::max_num_args_called > m->numDims();  
      bool direction_check = traits_t::direction >=  m->numDims(); 
      if(callable_check || direction_check) throw std::runtime_error("diffops setMesh: # of args in callables must be <= # of dims in mesh and direction must be < # of dims."); 
      
      using Evaluator = fornfdm::linops::internal::Evaluator<Derived>; 
      Evaluator eval(derived()); 

      // handle m_prod_before / m_prod_after logic against num args in callables
      if constexpr(traits_t::max_num_args_called == 0){
        // we can just store the compressed 1 Dimensional operator -> double kronecker product into correct dimension 
        m_prod_before = m->sizesMiddleProduct(0,traits_t::direction); 
        m_prod_after = m->sizesMiddleProduct(traits_t::direction+1, m->numDims());

        // resize + reserve the stencil 
        std::size_t nnz = Evaluator::nonZerosEstimate(axis); 
        m_stencil.reserve(nnz); 
        m_stencil.resize(axis_size, axis_size); 
        
        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          // use node selector 
          typename Evaluator::Row row(eval, m, m_prod_before * row_idx); 
          // copy the indices into m_stencil's inner indices ptr
          m_stencil.outerIndexPtr()[row_idx] = row.valuePtrOffset(); 
          std::copy_n(row.columnIndices().cbegin(), row.size(), m_stencil.innerIndexPtr() + row.valuePtrOffset());  
          row.mapToEigen(m_stencil.valuePtr() + row.valuePtrOffset()) = row.values(); 
        }
        // very last outer index needs set. 
        m_stencil.outerIndexPtr()[axis_size] = nnz; 
      }
      else if constexpr(traits_t::max_num_args_called <= traits_t::direction+1){
        // we can store the inflated 1st kronecker product, accounting for callables requiring coordinates. 
        m_prod_before = 1; 
        std::size_t product_before = m->sizesMiddleProduct(0,traits_t::direction); 
        m_prod_after = m->sizesMiddleProduct(traits_t::direction+1, m->numDims()); 

        // resize + reserve the stencil 
        std::size_t nnz = product_before * Evaluator::nonZerosEstimate(axis); 
        m_stencil.reserve(nnz); 
        m_stencil.resize(product_before*axis_size, product_before*axis_size); 

        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx < product_before*axis_size; ++row_idx)
        {
          // use node selector 
          typename Evaluator::Row row(eval, m, row_idx);
          // // copy the indices into m_stencil's inner indices ptr
          std::size_t inner_offset = (row.valuePtrOffset() * product_before)+(row.size()*(row_idx%product_before)); 
          std::size_t inset = row_idx%product_before; 
          std::transform(
            row.columnIndices().cbegin(), 
            std::next(row.columnIndices().cbegin(), row.size()), 
            m_stencil.innerIndexPtr() + inner_offset, 
            [&](std::size_t idx){ return inset + product_before*idx; }
          ); 
          m_stencil.outerIndexPtr()[row_idx] = inner_offset; 
          row.mapToEigen(m_stencil.valuePtr() + inner_offset) = row.values(); 
        }
        // very last outer index needs set. 
        m_stencil.outerIndexPtr()[product_before*axis_size] = nnz; 
      }
      else{ // direction+1 < max_num_args <= m->numDims()
        // we can store the inflated 1st kronecker product repeated n times 
        m_prod_before = 1; 
        std::size_t product_before = m->sizesMiddleProduct(0,traits_t::direction); 
        std::size_t num_repeats = m->sizesMiddleProduct(traits_t::direction+1, traits_t::max_num_args_called); 
        m_prod_after = m->sizesMiddleProduct(traits_t::max_num_args_called, m->numDims()); 

        std::size_t nnz = product_before * Evaluator::nonZerosEstimate(axis); 
        m_stencil.reserve(num_repeats * nnz); 
        m_stencil.resize(num_repeats*product_before*axis_size, num_repeats*product_before*axis_size); 

        // write node wise expressions into each row of stencil 
        for(std::size_t row_idx=0; row_idx < num_repeats*product_before*axis_size; ++row_idx)
        {
          typename Evaluator::Row row(eval, m, row_idx); 
          std::size_t inner_offset = (nnz)*(row_idx/(product_before*axis_size))+(row.valuePtrOffset() * product_before)+(row.size()*(row_idx%product_before)); 
          std::size_t inset = (product_before*axis_size)*(row_idx/(product_before*axis_size)) + (row_idx%product_before); 
          std::transform(
            row.columnIndices().cbegin(), 
            std::next(row.columnIndices().cbegin(), row.size()), 
            m_stencil.innerIndexPtr() + inner_offset, 
            [&](std::size_t idx){ return inset + product_before*idx; }
          ); 
          m_stencil.outerIndexPtr()[row_idx] = inner_offset; 
          row.mapToEigen(m_stencil.valuePtr() + inner_offset) = row.values();
        }
        // very last outer index needs set. 
        m_stencil.outerIndexPtr()[num_repeats*product_before*axis_size] = num_repeats * nnz; 
      } 
    }
}; 

} // end namespace linops
} // end namespace fornfdm 

#endif // PartialDerivBase.hpp  