// PartialDerivBase.hpp
//
// CRTP base class all Partial Derivatives and expressions will derive from.
// Expressions are only able to be constructed if they add, minus, scalar mult, etc...
// on the same "selected nodes" in space. (Nwise...Op ~ Node Wise binar/unar Operation)
// and are the same direction partial derivative.  
//
// JAF 4/17/2026 

#ifndef PARTIALDERIVBASE_H
#define PARTIALDERIVBASE_H

#include "EvaluatorBase.hpp"

namespace fdm{
namespace linops{

template<class Derived> class PartialDerivBase; 

namespace internal{ 

// Evaluator 
template<class Derived>
struct Evaluator<PartialDerivBase<Derived>> : public EvaluatorBase<PartialDerivBase<Derived>>
{
  Evaluator<Derived> m_derived_eval; 

  Evaluator(const PartialDerivBase<Derived>& xpr) : m_derived_eval(xpr.derived()){}

  template<std::size_t N>
  auto evaluateWeightsAndCoords(const fdm::Scalar* weights, std::size_t weights_per_order, const Coordinate<N>& coords)
  {
    return m_derived_eval.evaluateWeightsAndCoords(weights, weights_per_order, coords); 
  }
}; 


// linops traits
template<class Derived>
struct traits_impl<fdm::linops::PartialDerivBase<Derived>> : traits_impl<Derived>{}; 

} // end namespace internal 
} // end namespace internal 
} // end namespace fdm 

// Eigen traits 
namespace Eigen{
namespace internal{

// traits of PartialDerivBase is same as Derived
template<class Derived>
struct traits<fdm::linops::PartialDerivBase<Derived>> : public traits<Derived>{}; 

} // end namespace internal 
} // end namespace Eigen

namespace fdm{
  namespace linops{

template<class Derived>
class PartialDerivBase : public Eigen::SparseMatrixBase<Derived> // TODO inherit from compressed base? whats the difference?  
{
  public:
    // Type Defs --------------------- 
    typedef Eigen::SparseMatrixBase<Derived> Base;
    EIGEN_SPARSE_PUBLIC_INTERFACE(PartialDerivBase) 

    // Friends
    friend Eigen::internal::evaluator<PartialDerivBase>; 

  public: // TODO make private 
    // Member Data ----------------------------------------------
    // const Mesh* m_mesh_raw = nullptr; // This is unused until I need time dependent operators....... 
    std::weak_ptr<const Mesh> m_mesh_observed = {/*nullptr*/}; 
    fdm::CSRMatrix m_stencil = {}; 
    std::size_t m_prod_before = 1; 
    std::size_t m_prod_after = 1; 

  public:
    // Member Functions =====================================================
    // FDM Interface -------
    void setMesh(const std::shared_ptr<const Mesh>& m)
    {
      using traits_t = fdm::linops::internal::traits<Derived>; 
      const auto& axis = m->getAxis(traits_t::direction); 
      const std::size_t axis_size = m->sizeOfDim(traits_t::direction); 

      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      bool callable_check = traits_t::max_num_args_called > m->numDims();  
      bool direction_check = traits_t::direction >=  m->numDims(); 
      if(callable_check || direction_check) throw std::runtime_error("diffops setMesh: # of args in callables must be <= # of dims in mesh and direction must be < # of dims."); 
      
      // update state 
      this->m_mesh_observed = m; // does not hook! getMesh() returns nullptr on leafs ( they will never be calculated in this function)

      using Evaluator = fdm::linops::internal::Evaluator<Derived>; 
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
          typename Evaluator::Row row(eval, m.get(), m_prod_before * row_idx); 

          std::cout << "row: " << row_idx << " valuePtrOffset: " << row.valuePtrOffset() << std::endl; 

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

        std::cout << "setMesh: innerIndexPtr: " << m_stencil.innerIndexPtr() << std::endl;

        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx < product_before*axis_size; ++row_idx)
        {
          // use node selector 
          typename Evaluator::Row row(eval, m.get(), row_idx); 

          std::cout << "row: " << row_idx << " valuePtrOffset: " << (row.valuePtrOffset() * product_before)+(row.size()*(row_idx%product_before)) << std::endl; 

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
        m_stencil.outerIndexPtr()[product_before * axis_size] = nnz; 
      }
      else{ // direction+1 < max_num_args <= m->numDims()
        m_prod_before = 1; 
        m_prod_after = m->sizesMiddleProduct(traits_t::max_num_args_called, m->numDims()); 
        // this could be combined with the above implementation??????? 
        // TODO 
      }      
    }

    SharedConstMesh getMesh() const { return m_mesh_observed.lock(); }
    void setTime(double t){ derived().setTime(t); }
    double getTime() const { return derived().getTime(); }
    
    // Eigen Interface ------- 
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
}; 

  } // end namespace linops
} // end namespace fdm 

#endif // DiffOpBase.hpp  