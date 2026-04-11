// BlockDiagExpr.hpp
//
// take a sparse matrix and returns an Eigen expression
// representing kroenecker poroduct I(n) @ A
// where I(n) is identity with n rows   
//
// 2/7/2026 

#ifndef BLOCKDIAGEXPR_H
#define BLOCKDIAGEXPR_H

#include<Eigen/Sparse>

// Forward declarations ---------------------------------------------
namespace fdm{
  namespace utils{

template<typename ArgTpe>
class BlockDiag;

  }
}


// type traits =======================================================================
namespace Eigen {
namespace internal {
template<class ArgType>
struct traits<fdm::utils::BlockDiag<ArgType> > {
  typedef Eigen::Sparse StorageKind;
  typedef Eigen::MatrixXpr XprKind;
  typedef typename ArgType::StorageIndex StorageIndex;
  typedef typename ArgType::Scalar Scalar;
  enum {
    Flags = Eigen::RowMajor,
    RowsAtCompileTime = Eigen::Dynamic,
    ColsAtCompileTime = Eigen::Dynamic,
    MaxRowsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };
};
}  // namespace internal
}  // namespace Eigen

// expression class ======================================================================= 
namespace fdm{
  namespace utils{

template<class ArgType>
class BlockDiag : public Eigen::SparseMatrixBase< BlockDiag<ArgType> > {
  public:
    // typedefs 
    typedef typename Eigen::internal::ref_selector<BlockDiag>::type Nested;
    typedef Eigen::Index Index;
    typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
    
    // constructors 
    BlockDiag(const ArgType& arg_init, std::size_t repeats_init=1) : m_arg(arg_init), m_num_repeats(repeats_init) {
      // ArgType is RowMajor... 
      static_assert(Eigen::internal::traits<ArgType>::Flags & Eigen::RowMajorBit, "ArgType to BlockDiag must be RowMajor"); 
      if(repeats_init<1) throw std::invalid_argument("BlockDiag must be constructed with num_repeats >= 1."); 
    }
    
    // member functions 
    Index rows() const { return m_num_repeats * m_arg.rows(); }
    Index cols() const { return m_num_repeats * m_arg.cols(); }

    // member data 
    ArgTypeNested m_arg;
    std::size_t m_num_repeats; 
  
};

  } // end namespace utils 
} // end namespace fdm 

// the evaluator =======================================================================
namespace Eigen {
namespace internal {
template<typename ArgType>
struct evaluator< fdm::utils::BlockDiag<ArgType> > : evaluator_base< fdm::utils::BlockDiag<ArgType> > {

  // typedefs -------------------------------------------------- 
  typedef fdm::utils::BlockDiag<ArgType> XprType;
  typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
  typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // custom InnerIterator ----------------------------------
  struct InnerIterator{
    // Constructor ================================================================
    InnerIterator(const evaluator& eval, Index row_idx)
      : m_eval(eval), 
      m_row(row_idx), 
      m_col_offset(eval.m_xpr.m_arg.cols() * (row_idx / eval.m_xpr.m_arg.rows())), 
      m_wrapped_it(eval.m_argImpl, row_idx % eval.m_xpr.m_arg.rows())
    {};

    // Member Funcs ===================================================
    operator bool() const { return m_wrapped_it; }
    void operator++(){ ++m_wrapped_it; }
    Index row() const { return m_row; }
    Index col() const { return m_col_offset + m_wrapped_it.col(); }
    Index index() const { return m_col_offset + m_wrapped_it.col(); }
    Scalar value() const { return m_wrapped_it.value(); }
    // member data ------------------------------------------
    const evaluator& m_eval; 
    typename evaluator<ArgTypeNestedCleaned>::InnerIterator m_wrapped_it;
    Index m_row; 
    Index m_col_offset; 
  }; // end InnerIterator 

  // Constructors ======================================================== 
  evaluator(const XprType& xpr) 
    : m_argImpl(xpr.m_arg), m_xpr(xpr), m_num_repeats(xpr.m_num_repeats)
  {};
 
  // Member Functions ========================================================
  Index rows() const {return m_xpr.rows() * m_num_repeats; }; 
  Index cols() const {return m_xpr.cols() * m_num_repeats; }; 
  Index outerSize() const { return m_xpr.rows() * m_num_repeats; }
  Index innerSize() const { return m_xpr.cols() * m_num_repeats; }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate() * m_num_repeats; }

  // Flags ------------------------------------------------------
  enum { CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost, Flags = Eigen::RowMajor };
 
  // Member Data ------------------------------------------------------
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  const XprType& m_xpr;  
  Index m_num_repeats; 
};
}  // namespace internal
}  // namespace Eigen

// the entry point ======================================================================= 
namespace fdm{
  namespace utils{

template<class ArgType>
BlockDiag<ArgType> make_BlockDiag(const Eigen::SparseMatrixBase<ArgType>& arg, std::size_t num_repeats=1) {
  return BlockDiag<ArgType>(arg.derived(), num_repeats);
}

  } // end namespace utils 
} // end namespace fdm 

#endif // BlockDiagExpr.hpp 

