// FillRowExpr.hpp 
//
// Binary Expression f(A,B) that
// returns A's row iterator if it has non zero entries, 
// otherwise it returns B's row iterator
//
// JAF 3/27/2026 

#ifndef FILLROWSEXPR_H
#define FILLROWSEXPR_H 

namespace fdm{
namespace utils{ 

// Forward declarations ---------------------------------------------
template<typename ArgType01, typename ArgType02>
class FillRow;

} // end namespace utils 
} // end namespace fdm 

// type traits =======================================================================
namespace Eigen {
namespace internal {
template<typename ArgType01, typename ArgType02>
struct traits< fdm::utils::FillRow<ArgType01, ArgType02> > {
  typedef Eigen::Sparse StorageKind;
  typedef Eigen::MatrixXpr XprKind;
  typedef typename ArgType01::StorageIndex StorageIndex;
  typedef typename ArgType01::Scalar Scalar;
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

namespace fdm{
namespace utils { 

// expression class ======================================================================= 
template<typename ArgType01, typename ArgType02>
class FillRow : public Eigen::SparseMatrixBase< FillRow<ArgType01,ArgType02> > {
  public:
    // typedefs 
    typedef typename Eigen::internal::ref_selector<FillRow>::type Nested;
    typedef Eigen::Index Index;
    typedef typename Eigen::internal::ref_selector<ArgType01>::type ArgTypeNested01;
    typedef typename Eigen::internal::ref_selector<ArgType02>::type ArgTypeNested02;
    
    // constructors 
    FillRow(const ArgType01& A,const ArgType02& B) : m_arg01(A), m_arg02(B) {
      // Both ArgType01 (02) is RowMajor... 
      static_assert(Eigen::internal::traits<ArgType01>::Flags & Eigen::internal::traits<ArgType02>::Flags & Eigen::RowMajorBit, "ArgType to FillRow must be RowMajor"); 
      if(A.rows()!=B.rows() || A.cols()!=B.cols()) throw std::runtime_error("Error in FillRow constructo rows/cols don't match!"); 
    }; 
    
    // member functions 
    Index rows() const { return m_arg01.rows(); }
    Index cols() const { return m_arg01.cols(); }

    // member data 
    ArgTypeNested01 m_arg01;
    ArgTypeNested02 m_arg02;
};

} // end namespace utils 
} // end namespace fdm 

// the evaluator =======================================================================
namespace Eigen {
namespace internal {
template<typename ArgType01, typename ArgType02>
struct evaluator< fdm::utils::FillRow<ArgType01,ArgType02> > : evaluator_base< fdm::utils::FillRow<ArgType01,ArgType02> > {

  // typedefs -------------------------------------------------- 
  typedef fdm::utils::FillRow<ArgType01,ArgType02> XprType;
  typedef typename nested_eval<ArgType01, XprType::ColsAtCompileTime>::type ArgTypeNested01;
  typedef typename nested_eval<ArgType02, XprType::ColsAtCompileTime>::type ArgTypeNested02;
  // using ArgTypeNested = ArgTypeNested01; // does eigen need to have the ArgTypeNested type defined?  
  typedef typename remove_all<ArgTypeNested01>::type ArgTypeNestedCleaned01;
  typedef typename remove_all<ArgTypeNested02>::type ArgTypeNestedCleaned02;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // custom InnerIterator ----------------------------------
  struct InnerIterator{
    // Constructor ================================================================
    InnerIterator(const evaluator& eval, Index row_idx)
      : m_eval(eval),
      m_row(row_idx), 
      m_wrapped_it01(eval.m_argImp01, row_idx),
      m_wrapped_it02(eval.m_argImp02, row_idx)
    {
      m_use01 = (m_wrapped_it01); // AFTER m_wrapped_it01 is constructed ... 
    };

    // Member Funcs ===================================================
    operator bool() const { return m_use01 ? static_cast<bool>(m_wrapped_it01) : static_cast<bool>(m_wrapped_it02); }
    void operator++(){ 
      if(m_use01){ ++m_wrapped_it01; } else { ++m_wrapped_it02; }; 
    }
    Index row() const { return m_row; }
    Index col() const { return m_use01 ? m_wrapped_it01.col() : m_wrapped_it02.col(); }
    Index index() const { return m_use01 ? m_wrapped_it01.index() : m_wrapped_it02.index(); }
    Scalar value() const { return m_use01 ? m_wrapped_it01.value() : m_wrapped_it02.value(); }
    // member data ------------------------------------------
    const evaluator& m_eval; 
    bool m_use01; 
    typename evaluator<ArgTypeNestedCleaned01>::InnerIterator m_wrapped_it01;
    typename evaluator<ArgTypeNestedCleaned02>::InnerIterator m_wrapped_it02; 
    Index m_row; 

  }; // end InnerIterator 

  // Constructors ======================================================== 
  evaluator(const XprType& xpr) 
    : m_argImp01(xpr.m_arg01), m_argImp02(xpr.m_arg02), m_xpr(xpr)
  {}
 
  // Member Functions ========================================================
  Index rows() const {return m_xpr.rows(); } 
  Index cols() const {return m_xpr.cols(); }
  Index outerSize() const { return m_xpr.rows(); }
  Index innerSize() const { return m_xpr.cols(); }
  Index nonZerosEstimate() const { return m_xpr.nonZerosEstimate(); }

  // Flags ------------------------------------------------------
  enum { 
    CoeffReadCost = evaluator<ArgTypeNestedCleaned01>::CoeffReadCost, 
    Flags = Eigen::RowMajor
  };
 
  // Member Data ------------------------------------------------------
  evaluator<ArgTypeNestedCleaned01> m_argImp01;
  evaluator<ArgTypeNestedCleaned02> m_argImp02;
  const XprType& m_xpr;  
};
}  // namespace internal
}  // namespace Eigen

namespace fdm{
namespace utils{ 

// the entry point ======================================================================= 
template<class ArgType01, class ArgType02>
FillRow<ArgType01,ArgType02> make_FillRow(const Eigen::SparseMatrixBase<ArgType01>& A,const Eigen::SparseMatrixBase<ArgType02>& B) {
  return FillRow<ArgType01, ArgType02>(A.derived(), B.derived());
}

template<class T, class U>
auto make_FillRow_fold(T&& A, U&& B){
  return make_FillRow(std::forward<T>(A), std::forward<U>(B)); 
}

template<class T, class U, typename... Args>
auto make_FillRow_fold(T&& A, U&& B, Args&&... rest){
  auto combined = make_FillRow(std::forward<T>(A), std::forward<U>(B)); 
  
  if constexpr(sizeof...(Args)==0){
    return combined; 
  }
  else{
    return make_FillRow_fold(std::move(combined), std::forward<Args>(rest)...); 
  }
}

} // end namespac utils 
} // end namespace fdm

#endif 