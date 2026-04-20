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

#include<FiniteDifference/Utilities/FornbergArrayCalc.hpp> 
#include "NodeSelector.hpp"
#include "CoordinateSelector.hpp" 

namespace fdm{
namespace linops{
template<class Derived>
class PartialDerivBase; 
}
namespace internal{

// FDM traits
template<class Derived>
struct traits_impl<fdm::linops::PartialDerivBase<Derived>> : traits_impl<Derived>{}; 

} // end namespace internal 
} // end namespace fdm 

// Eigen traits 
namespace Eigen{
  namespace internal{

// all Partial Derivatives + Expressions will act as the same Eigen traits 
// wait eigen has a nestbyref bit inside of traits we can resuse that.......  
template<class Derived>
struct traits<fdm::linops::PartialDerivBase<Derived>>
{
  typedef fdm::Scalar Scalar;
  typedef Eigen::Index StorageIndex;
  typedef Sparse StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = Eigen::RowMajor | NestByRefBit /* | no assignment LvalueBit  */ /* | not CompressedAccessBit*/ ,
    SupportedAccessPatterns = OuterRandomAccessPattern
  };
}; 

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
    // typedef linops::internal::NodeSelector<Derived> NodeSelector; 
    template<std::size_t N> using NodeSelector = fdm::linops::internal::NodeSelector<N>; 
    // typedef linops::internal::RowEvaluator<Derived> RowEvaluator; // TODO this can go outside the class? 
    typedef CSRMatrix::Index Index; 

    // Friends
    friend Eigen::internal::evaluator<PartialDerivBase>; 
    // friend fdm::internal::evaluator<PartialDerivBase>; TODO name the fdm linops evaluator 

  public: // TODO make private 
    // Member Data ----------------------------------------------
    // const Mesh* m_mesh_raw = nullptr; // This is unused until I need time dependent operators....... 
    std::weak_ptr<const Mesh> m_mesh_observed = {/*nullptr*/}; 
    fdm::CSRMatrix m_stencil = {}; 
    std::size_t m_prod_before = 1; 
    std::size_t m_prod_after = 1; 

  public:
    // Member Functions =====================================================
    template<std::size_t numNodesMax, std::size_t numCoordsMax=0>
    auto evaluateWeightsAndCoords(
      const std::array<double, numNodesMax>& weights, 
      std::size_t weights_per_order, 
      const std::array<double,numCoordsMax>& coords={}) const 
    { return derived().evaluateWeightsAndCoords(weights, weights_per_order, coords); }

    // FDM Interface -------
    void setMesh(const std::shared_ptr<const Mesh>& m)
    {
      constexpr std::size_t max_num_args_called = fdm::internal::traits<Derived>::max_num_args_called; 
      constexpr int direction = fdm::internal::traits<Derived>::direction;
      constexpr std::size_t max_order = fdm::internal::traits<Derived>::maxOrder; 
      const auto& axis = m->getAxis(direction); 
      const std::size_t axis_size = m->sizeOfDim(direction); 

      // if number of args in callable c(x,y,z) is > meshes # of dims throw 
      bool callable_check = max_num_args_called > m->numDims();  
      bool direction_check = direction >=  m->numDims(); 
      if(callable_check || direction_check) throw std::runtime_error("diffops setMesh: # of args in callables must be <= # of dims in mesh and direction must be < # of dims."); 
      
      // update state 
      this->m_mesh_observed = m; // does not hook! getMesh() returns nullptr on leafs ( they will never be calculated in this function)
      // all pieces in expression tree are looking at the same Mesh. ready to use selectors + evaluators  

      // handle m_prod_before / m_prod_after logic against num args in callables
      if constexpr(max_num_args_called == 0){
        // we can just store the compressed 1 Dimensional operator -> double kronecker product into correct dimension 
        m_prod_before = m->sizesMiddleProduct(0,direction); 
        m_prod_after = m->sizesMiddleProduct(direction+1, m->numDims());

        // resize + reserve the stencil 
        std::size_t nnz = fdm::linops::internal::NodeSelector< max_order+1 >::sumNodesPerRow(m->getAxis(direction)); 
        m_stencil.reserve(nnz); 
        m_stencil.resize(axis_size, axis_size); 
        
        // write node wise expressions into each row stencil 
        for(std::size_t row_idx=0; row_idx<axis_size; ++row_idx)
        {
          // use node selector 
          const NodeSelector<max_order+1> nodes(axis, row_idx);

          // if max_num_args_called > 0 this would use coordinates too. // const CoordinateSelector coord_selector(m, row_idx);  

          // Fornberg algorithm on the stack. 
          fdm::utils::FornArrayCalc<nodes.numNodesMax, max_order> weight_calc; 
          weight_calc.calculate(nodes.x_bar, nodes.nodeValues.cbegin(), std::next(nodes.nodeValues.cbegin(), nodes.numNodesUsed)); 

          // copy the indices into m_stencil's inner indices ptr 
          m_stencil.outerIndexPtr()[row_idx] = nodes.nonZerosOffset; 
          std::copy_n(nodes.nodeIndices.begin(), nodes.numNodesUsed, m_stencil.innerIndexPtr()+nodes.nonZerosOffset); 

          // copy the derived's expression into m_stencil non zeros. Eigen::Map handles SIMD :-)  
          using Mapped = Eigen::Map<Eigen::Matrix<fdm::Scalar, 1, Eigen::Dynamic>>;  
          Mapped( m_stencil.valuePtr()+nodes.nonZerosOffset, nodes.numNodesUsed ) = this->evaluateWeightsAndCoords(weight_calc.getArray(),nodes.numNodesUsed);          
        }
        // veeeerrry last outer index needs set. 
        m_stencil.outerIndexPtr()[axis_size] = nnz; 
      }
      else if constexpr(max_num_args_called <= direction){
        // we can store the inflated 1st kronecker product, accounting for callables requiring coordinates. 
        m_prod_before = 1; 
        m_prod_after = m->sizesMiddleProduct(direction+1, m->numDims()); 
      }
      else{ // max_num_args <= m->numDims()
        m_prod_before = 1; 
        m_prod_after = m->sizesMiddleProduct(max_num_args_called, m->numDims()); 
        // this could be combined with the above implementation??????? 
      }      
    }

    SharedConstMesh getMesh() const { try{ return m_mesh_observed.lock(); } catch(...){ return nullptr; } }
    void setTime(double t){ derived().setTime(t); }
    void getTime() const { return derived().getTime(); }
    
    // Eigen Interface ------- 
    Index rows() const { return m_prod_before * m_prod_after * m_stencil.rows(); }
    Index cols() const { return m_prod_before * m_prod_after * m_stencil.cols(); }
    Index nonZerosEstimate() const {return m_prod_before * m_prod_after * m_stencil.nonZeros(); }

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