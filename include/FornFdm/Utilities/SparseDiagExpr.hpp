// SparseDiagExpr.hpp
//
// takes a Eigen::SparseMatrix<double,Eigen::RowMajor> matrix of rows==1 and makes it into a diagonal matrix 
// tightly based on eigen docs for make circulant
// https://libeigen.gitlab.io/eigen/docs-5.0/TopicNewExpressionType.html
// 
// JAF 1/2/2026 

#ifndef FORNFDM_UTILS_SPARSEDIAGEXPR_H
#define FORNFDM_UTILS_SPARSEDIAGEXPR_H

#include<cstdint>
#include<Eigen/Core>
#include<Eigen/Sparse>

namespace fornfdm{
  namespace utils{

// Enum + Forward declarations ---------------------------------------------
enum class SparseDiagPattern { REPEAT, CYCLE }; 
template<typename ArgTpe, SparseDiagPattern P>
class SparseDiag; 

  } // end namespace utils 
} // end namespace fornfdm 

// type traits =======================================================================
namespace Eigen {
namespace internal {
template<class ArgType, fornfdm::utils::SparseDiagPattern P>
struct traits<fornfdm::utils::SparseDiag<ArgType, P> > {
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
namespace fornfdm{
  namespace utils{

template<class ArgType, SparseDiagPattern P = SparseDiagPattern::REPEAT>
class SparseDiag : public Eigen::SparseMatrixBase< SparseDiag<ArgType,P> > {
  public:
    // typedefs 
    typedef typename Eigen::internal::ref_selector<SparseDiag>::type Nested;
    typedef Eigen::Index Index;
    typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
    
    // constructors 
    SparseDiag(const ArgType& arg_init, std::size_t repeats_init=1) : m_arg(arg_init), m_num_repeats(repeats_init) {
      // static_assert( std::is_same<ArgType, Matrix>::value, "Invalid arg type. SparseDiag must be constructed from Matrix."); 
      if(arg_init.rows()!=1 && arg_init.cols()!=1) throw std::invalid_argument("SparseDiag must be constructed from expr with .rows()==1 or .cols()==1"); 
      if(repeats_init<1) throw std::invalid_argument("SparseDiag must be constructed with num_repeats >= 1."); 
      m_is_horizontal = arg_init.rows()==1; 
      m_size = std::max(arg_init.rows(), arg_init.cols()); 
    }
    
    // member functions 
    Index rows() const { return m_num_repeats * m_size; }
    Index cols() const { return m_num_repeats * m_size; }

    // member data 
    ArgTypeNested m_arg;
    std::size_t m_num_repeats; 
    std::size_t m_size; 
    bool m_is_horizontal; 
  
};

  } // end namespace utils 
} // end namespace fornfdm 


// the evaluator =======================================================================
namespace Eigen {
namespace internal {
template<typename ArgType, fornfdm::utils::SparseDiagPattern P>
struct evaluator< fornfdm::utils::SparseDiag<ArgType,P> > : evaluator_base< fornfdm::utils::SparseDiag<ArgType,P> > {

  // typedefs -------------------------------------------------- 
  typedef fornfdm::utils::SparseDiag<ArgType, P> XprType;
  typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
  typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // custom InnerIterator ----------------------------------
  struct InnerIterator{
    // Constructor ================================================================
    InnerIterator(const evaluator& eval, Index row_idx)
      : m_eval(eval), m_row(row_idx)
    {
      // iterator is only valid if row_idx is in bounds and underlying row has nonzero entry 
      m_valid = (m_row < m_eval.rows()) && \
                      (value() != Scalar(0.0)); 
    };
    // Member Funcs ===================================================
    operator bool() const { return m_valid; }; 
    void operator++(){ m_valid=false; }; 
    Index row() const { return m_row; }; 
    Index col() const { return m_row; };
    Index index() const { return m_row; };
    Scalar value() const { 
      if constexpr(P == fornfdm::utils::SparseDiagPattern::REPEAT){
        if(m_eval.m_is_horizontal) return m_eval.m_argImpl.coeff(0, m_row / m_eval.m_num_repeats); 
        else return m_eval.m_argImpl.coeff(m_row / m_eval.m_num_repeats, 0);
      }
      if constexpr(P == fornfdm::utils::SparseDiagPattern::CYCLE){
        if(m_eval.m_is_horizontal) return m_eval.m_argImpl.coeff(0, m_row % (m_eval.rows()/m_eval.m_num_repeats));  
        else return m_eval.m_argImpl.coeff(m_row % (m_eval.rows()/m_eval.m_num_repeats), 0); 
      }
      else{
        if(m_eval.m_is_horizontal) return m_eval.m_argImpl.coeff(0, m_row / m_eval.m_num_repeats); 
        else return m_eval.m_argImpl.coeff(m_row / m_eval.m_num_repeats, 0);
      }
    };
    // member data ------------------------------------------
    const evaluator& m_eval; 
    Index m_row;    
    bool m_valid; 
  }; // end InnerIterator 

  // Constructors ======================================================== 
  evaluator(const XprType& xpr) 
    : m_argImpl(xpr.m_arg), m_rows(xpr.rows()), m_cols(xpr.cols()), 
    m_num_repeats(xpr.m_num_repeats), m_is_horizontal(xpr.m_is_horizontal) 
  {};
 
  // Member Functions ========================================================
  Index rows() const {return m_rows; }; 
  Index cols() const {return m_cols; }; 
  Index outerSize() const { return m_rows; }
  Index innerSize() const { return m_rows; }
  Index nonZerosEstimate() const { return m_rows; }

  // Flags ------------------------------------------------------
  enum { CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost, Flags = Eigen::RowMajor };
 
  // Member Data ------------------------------------------------------
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  Index m_rows; 
  Index m_cols; 
  Index m_num_repeats; 
  Index m_size; 
  bool m_is_horizontal; 
};
}  // namespace internal
}  // namespace Eigen

// the entry point ======================================================================= 
namespace fornfdm{
  namespace utils{

template<SparseDiagPattern P = SparseDiagPattern::REPEAT, class ArgType>
SparseDiag<ArgType,P> make_SparseDiag(const Eigen::SparseMatrixBase<ArgType>& arg, std::size_t num_repeats=1) {
  return SparseDiag<ArgType, P>(arg.derived(), num_repeats);
}

template<SparseDiagPattern P = SparseDiagPattern::REPEAT, class ArgType>
SparseDiag<ArgType,P> make_SparseDiag(const Eigen::MatrixBase<ArgType>& arg, std::size_t num_repeats=1) {
  return SparseDiag<ArgType, P>(arg.derived(), num_repeats);
}

  } // end namespace utils 
} // end namespace fornfdm 

#endif // SparseDiagExpr.hpp