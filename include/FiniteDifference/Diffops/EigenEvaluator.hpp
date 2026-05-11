// EigenEvaluator.hpp 
//
// Handles reading the CSRMatrix inside of a PartialDerivBase<> 
// into a kronecker product eigen expression  
//
// JAF 4/17/2026 

#ifndef FDM_DIFFOPS_EIGENEVALUATOR_H
#define FDM_DIFFOPS_EIGENEVALUATOR_H

#include<FiniteDifference/Utilities/BlockDiagExpr.hpp>
#include<FiniteDifference/Utilities/HighDimExpr.hpp>

namespace Eigen {
namespace internal {

// // implementation depends on if max_num_args_called == 0, <= maxOrder, or > maxOrder 
// template<class Derived, std::size_t max_num_args_called>
// struct EigenEvaluator_impl; 

// reads row iterators from double kronecker products I(n) D @ I(m)
template<class Derived>
struct evaluator<fdm::linops::PartialDerivBase<Derived>> : public evaluator_base<fdm::linops::PartialDerivBase<Derived>>
{
  typedef Derived XprType;

  typedef typename nested_eval< fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>> , XprType::ColsAtCompileTime >::type ArgTypeNested;
  typedef typename remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; // do we need theeeesee? 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 
  enum { CoeffReadCost = evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>>::CoeffReadCost, Flags = Eigen::RowMajorBit };

  struct InnerIterator : public evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>>::InnerIterator
  {
    InnerIterator(const evaluator& eval, Index row_idx)
      : evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>>::InnerIterator(eval.m_stencilKroneckerImpl, row_idx)
    {}
    InnerIterator(const Eigen::internal::evaluator<Derived>& eval, Index row_idx)
      : evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>>::InnerIterator(eval.m_stencilKroneckerImpl, row_idx)
    {}
  }; 

  evaluator(const fdm::linops::PartialDerivBase<Derived>& xpr_b)
    : m_xpr(xpr_b.derived()), 
    m_stencilKronecker(fdm::utils::make_BlockDiag(fdm::utils::make_HighDim(m_xpr.m_stencil, m_xpr.m_prod_before), m_xpr.m_prod_after)), 
    m_stencilKroneckerImpl(m_stencilKronecker)
  {}

  evaluator(const Derived& xpr)
    : m_xpr(xpr), 
    m_stencilKronecker(fdm::utils::make_BlockDiag(fdm::utils::make_HighDim(m_xpr.m_stencil, m_xpr.m_prod_before), m_xpr.m_prod_after)), 
    m_stencilKroneckerImpl(m_stencilKronecker)
  {}

  // Copy 
  evaluator(const evaluator& other)
    : m_xpr(other.m_xpr),
      m_stencilKronecker(fdm::utils::make_BlockDiag(fdm::utils::make_HighDim(m_xpr.m_stencil, m_xpr.m_prod_before), m_xpr.m_prod_after)),
      m_stencilKroneckerImpl(m_stencilKronecker) // Re-bind to local member!
  {}

  Index rows() const {return m_xpr.rows(); }
  Index cols() const {return m_xpr.cols(); }
  Index outerSize() const { return m_xpr.rows(); }
  Index innerSize() const { return m_xpr.cols(); }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // don't pay attention to NestByRefBit at all. m_stencil should never be copied. 
  const XprType& m_xpr; 
  const fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>> m_stencilKronecker; 
  evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>> m_stencilKroneckerImpl; 
}; 

// reads row iterators from single kronecker products D @ I(m)

}  // namespace internal
}  // namespace Eigen

#endif // EigenEvaluator.hpp

// template<class Derived>
// struct diffop_evaluator_base
// {
//   // Type Defs ----------------------------------
//   typedef diffop_evaluator_base<Derived> Base;
//   typedef diffop_evaluator_base<Derived>::DerivedAsNestedType NestedType;   
//   // Member Functions ------------------------------------------------------ 
//   static void reseatMap(Derived& eval, )
//     void reseatMap(const Calculator& c, const int& numNodes)
//   {
//     new (&m_vals) = Eigen::::Map< Eigen::VectorXd>( c.getArray().data() + order * numNodes, numNodes);     
//   }
//   // Nested Struct -------------------------------------------------- 
//   struct DerivedAsNestedType 
//   {
    
//   }

//   // CRTP casts 
//   inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
//   inline Derived& derived() { return *static_cast<Derived*>(this); }
//   // inline Derived& const_cast_derived() const /* not necessary? */  
//   // { return *static_cast<Derived*>(const_cast<SparseMatrixBase*>(this)); }

// } 