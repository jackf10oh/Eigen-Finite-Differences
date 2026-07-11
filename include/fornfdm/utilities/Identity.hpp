// Identity.hpp
//
// Row major expressions representing m x n identity
// 
// JAF 1/2/2026 

#ifndef FORNFDM_UTILS_IDENTITy_H
#define FORNFDM_UTILS_IDENTITy_H

#include<cstdint>
#include<Eigen/Sparse>
#include "../types.hpp"

// Enum + Forward declarations ---------------------------------------------
namespace fornfdm{
  namespace utils{
    class Identity; 
  }
}

// type traits =======================================================================
namespace Eigen {
namespace internal {
template<>
struct traits<fornfdm::utils::Identity> {
  typedef Eigen::Sparse StorageKind;
  typedef Eigen::MatrixXpr XprKind;
  typedef typename Eigen::Index StorageIndex;
  typedef fornfdm::Scalar Scalar;
  enum {
    Flags = Eigen::RowMajorBit,
    RowsAtCompileTime = Eigen::Dynamic,
    ColsAtCompileTime = Eigen::Dynamic,
    MaxRowsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };
};
}  // namespace internal
}  // namespace Eigen

// expression class ======================================================================= 
namespace fornfdm{
  namespace utils{

class Identity : public Eigen::SparseMatrixBase< Identity > {
  public:
    // typedefs 
    // typedefs 
    typedef typename Eigen::internal::ref_selector<fornfdm::utils::Identity>::type Nested;
    typedef Eigen::Index Index;
    
    // constructors 
    Identity(std::size_t rows, std::size_t cols) : m_rows(rows), m_cols(cols){}
    
    // member functions 
    Index rows() const { return m_rows; }
    Index cols() const { return m_cols; }
    void resize(std::size_t m, std::size_t n){m_rows=m; m_cols=n;}

    // member data 
    std::size_t m_rows; 
    std::size_t m_cols;   
};

  } // end namespace utils 
} // end namespace fornfdm 

// the evaluator =======================================================================
namespace Eigen {
  namespace internal {

template<>
struct evaluator< fornfdm::utils::Identity > : evaluator_base< fornfdm::utils::Identity > {

  // typedefs -------------------------------------------------- 
  typedef fornfdm::utils::Identity XprType;
  // typedef typename nested_eval<XprType::ColsAtCompileTime>::type ArgTypeNested;
  // typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // custom InnerIterator ----------------------------------
  struct InnerIterator{
    // Constructor ================================================================
    InnerIterator(const evaluator& eval, Index row_idx)
      : m_valid(row_idx < eval.cols()), m_row(row_idx)
    {};
    // Member Funcs ===================================================
    operator bool() const { return m_valid; } 
    void operator++(){ m_valid=false; }
    Index row() const { return m_row; }
    Index col() const { return m_row; }
    Index index() const { return m_row; }
    Scalar value() const { return 1.0; }

    // member data ------------------------------------------
    Index m_row;
    bool m_valid; 
  }; // end InnerIterator 

  // Constructors ======================================================== 
  evaluator(const XprType& xpr) 
    : m_rows(xpr.rows()), m_cols(xpr.cols()) 
  {};
 
  // Member Functions ========================================================
  Index rows() const {return m_rows; }; 
  Index cols() const {return m_cols; }; 
  Index outerSize() const { return m_rows; }
  Index innerSize() const { return m_cols; }
  Index nonZerosEstimate() const { return std::min(m_rows,m_cols); }

  // Flags ------------------------------------------------------
  enum { CoeffReadCost = 0, Flags = Eigen::RowMajor };
 
  // Member Data ------------------------------------------------------
  std::size_t m_rows; 
  std::size_t m_cols; 
};

  }  // namespace internal
}  // namespace Eigen

// the entry point ======================================================================= 
namespace fornfdm{
  namespace utils{

Identity make_Identity(std::size_t m, std::size_t n) {
  return Identity(m,n); 

    }
  }
}

#endif // Identity.hpp