#ifndef FOOTESTCLASS_H
#define FOOTESTCLASS_H

// #include "FiniteDifference/LinOps/LinOpTraits.hpp"
#include "FiniteDifference/Mesh.hpp"

struct Foo; 

namespace Eigen{
using FdmMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>; 

namespace internal{

template<>
struct traits<Foo> : public traits<FdmMatrix>{};

} // namespace internal
}// namespace Eigen 

class Foo : public Eigen::SparseCompressedBase<Foo> 
{
  private:
    // Type Defs ---------------------------------- 
  public:
    typedef Eigen::SparseCompressedBase<Foo> Base;
    using Base::convert_index;
    typedef typename Eigen::internal::traits<Foo>::Scalar Scalar; /*!< \brief Numeric type, e.g. float, double, int or std::complex<float>. */ \
    typedef typename Eigen::NumTraits<Scalar>::Real RealScalar; /*!< \brief The underlying numeric type for composed scalar types. \details In cases where Scalar is e.g. std::complex<T>, T were corresponding to RealScalar. */ \
    typedef typename Base::CoeffReturnType CoeffReturnType; /*!< \brief The return type for coefficient access. \details Depending on whether the object allows direct coefficient access (e.g. for a MatrixXd), this type is either 'const Scalar&' or simply 'Scalar' for objects that do not allow direct coefficient access. */ \
    typedef typename Eigen::internal::ref_selector<Foo>::type Nested; 
    typedef typename Eigen::internal::traits<Foo>::StorageKind StorageKind; 
    typedef typename Eigen::internal::traits<Foo>::StorageIndex StorageIndex; 
    enum CompileTimeTraits 
        { RowsAtCompileTime = Eigen::internal::traits<Foo>::RowsAtCompileTime, 
          ColsAtCompileTime = Eigen::internal::traits<Foo>::ColsAtCompileTime, 
          Flags = Eigen::internal::traits<Foo>::Flags, 
          SizeAtCompileTime = Base::SizeAtCompileTime, 
          MaxSizeAtCompileTime = Base::MaxSizeAtCompileTime, 
          IsVectorAtCompileTime = Base::IsVectorAtCompileTime }; 
    using Base::derived; 
    using Base::const_cast_derived;

    friend Eigen::internal::evaluator<Foo>; 
    // friend class SparseVector<_Scalar,0,_StorageIndex>;
    // template<typename, typename, typename, typename, typename> friend struct internal::Assignment; // no assignment allowed 

    // Constructors + Destructor ------------------------------ 
    Foo()=default; 
    Foo(const Foo& other) {
      std::cout << "Foo copied\n";
    }
    ~Foo()=default; 

    // Eigen Interface ---------------------------------- 
    auto rows() const { return m_mat.rows(); }
    auto cols() const { return m_mat.cols(); }
    void resize(Eigen::Index i, Eigen::Index j){ m_mat.resize(i,j); }
    inline const Scalar* valuePtr() const { return m_mat.valuePtr(); }
    inline Scalar* valuePtr() { return m_mat.valuePtr(); }
    inline const StorageIndex* innerIndexPtr() const { return m_mat.innerIndexPtr(); }
    inline StorageIndex* innerIndexPtr() { return m_mat.innerIndexPtr(); }
    inline const StorageIndex* outerIndexPtr() const { return m_mat.outerIndexPtr(); }
    inline StorageIndex* outerIndexPtr() { return m_mat.outerIndexPtr(); }
    inline const StorageIndex* innerNonZeroPtr() const { return m_mat.innerNonZeroPtr(); }
    inline StorageIndex* innerNonZeroPtr() { return m_mat.innerNonZeroPtr(); }

    // Deleted ----------------- 
    template<typename OtherDerived>
    inline auto& assign(const OtherDerived& other)=delete; 

    template<typename OtherDerived>
    inline void assignGeneric(const OtherDerived& other)=delete; 

    // FDM interface --------------------- 
    void setMesh(const std::shared_ptr<const fdm::Mesh>& m){ std::cout <<"Foo SetMesh!" << std::endl;  m_mesh_ptr=m; }
    auto getMesh()const{return m_mesh_ptr.lock(); }

    struct is_timedep_tag{}; // put this in derived classes if they override setTime() in a meaniful way. i.e. update state of linear operator
    void setTime(double t){ m_time = t; }
    double getTime()const{ return m_time; }

  private:
    // member data -----------------------------------
    Eigen::FdmMatrix m_mat; 
    std::weak_ptr<const fdm::Mesh> m_mesh_ptr; 
    double m_time; 

}; 

namespace Eigen{
namespace internal{

template<>
struct evaluator< Foo >
  : evaluator< FdmMatrix >
{
  typedef evaluator< FdmMatrix  > Base;
  evaluator() : Base() {}
  explicit evaluator(const Foo &f) : Base(f.m_mat) {}
};

}
}

#endif 