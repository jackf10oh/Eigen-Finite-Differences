// EigenEvaluator.hpp 
//
// Handles reading the CSRMatrix inside of a PartialDerivBase<> 
// into a kronecker product eigen expression  
//
// JAF 4/17/2026 

#ifndef DIFFOPSEIGENEVALUATOR_H
#define DIFFOPSEIGENEVALUATOR_H 

#include<FiniteDifference/Utilities/BlockDiagExpr.hpp>
#include<FiniteDifference/Utilities/HighDimExpr.hpp>

namespace fdm{
  namespace linops{
    namespace internal{

    } // endn namespace internal 
  } // end namespace linops 
} // end namespace fdm 

namespace Eigen {
namespace internal {
template<class Derived, std::size_t max_num_args_called>
struct EigenEvaluator_impl; 

template<class Derived>
struct evaluator<fdm::linops::PartialDerivBase<Derived>>: EigenEvaluator_impl< fdm::linops::PartialDerivBase<Derived>, fdm::internal::traits<Derived>::max_num_args_called >{}; 

// reads row iterators from double kronecker products I(n) D @ I(m) no need to pack coords
template<class Derived, std::size_t max_num_args_called, int direction>
struct EigenEvaluator_impl< fdm::linops::PartialDerivBase<Derived>,0>
{
  typedef typename Derived XprType;
  typedef typename nested_eval< fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>> , XprType::ColsAtCompileTime >::type ArgTypeNested;
  typedef typename remove_all< ArgTypeNested >::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType; // do we need theeeesee? 
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // Custom InnerIterator -----------------------------------------
  struct InnerIterator : public Eigen::internal::evaluator< fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>> >::InnerIterator
  {
    InnerIterator(const EigenEvaluator_impl& eval, Index row_idx) 
      : evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>>::InnerIterator(eval.m_stencilKronecker, row_idx)
    {} 
  };

  // Constructors ======================================================== 
  evaluator(const PartialDerivBase<Derived>& xpr) 
    : m_xpr(xpr.derived()), 
    m_stencilKronecker(fdm::utils::make_BlockDiag(fdm::utils::make_HighDim(m_xpr.m_stencil),m_xpr.m_prod_before)m_xpr.m_prod_after)
    // m_stencilKroneckerImpl(m_stencilKronecker)
  {};
  
  Index rows() const {return m_xpr.rows(); }; 
  Index cols() const {return m_xpr.cols(); }; 
  Index outerSize() const { return m_xpr.rows(); }
  Index innerSize() const { return m_xpr.cols(); }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Flags ------------------------------- 
  enum { CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost, Flags = Eigen::RowMajor };

  // Member Data ------------------ 
  const XprType& m_xpr; 
  const fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>> m_stencilKronecker; 
  // evaluator<fdm::utils::BlockDiag<fdm::utils::HighDim<fdm::CSRMatrix>>> m_stencilKroneckerImpl; // TODO unnecessary? 
}; 

// reads row iterators from single kronecker products D @ I(m)

}  // namespace internal
}  // namespace Eigen









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

#endif // Evaluator.hpp 