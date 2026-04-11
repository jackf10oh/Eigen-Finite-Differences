// LinOpTraits.hpp
//
// List of TMP functions for linops 
//
// JAF 12/7/2025

#ifndef LINOPTRAITS_H
#define LINOPTRAITS_H

#include<type_traits>
#include<Eigen/SparseCore>
#include<Eigen/Sparse> 
  
// Forward Declarations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace fdm{

class Vector1D; 
class VectorXD; 

class Mesh1D; 
using SharedConstMesh1D = std::shared_ptr<const Mesh1D>;
using WeakConstMesh1D = std::weak_ptr<const Mesh1D>;

class MeshXD; 
using SharedConstMeshXD = std::shared_ptr<const MeshXD>;
using WeakConstMeshXD = std::weak_ptr<const MeshXD>;

// conveniences - - - - - - - - - - - - -
using StridedRef = typename Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
using StrideView =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride>;
using Matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>; 

} // end namespace fdm 

namespace fdm{
  namespace linops{

template<typename Derived>
class LinOpMixIn;

template<typename Derived>
class LinOpBase1D;

template<typename Derived>
class LinOpBaseXD;

template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr;  

// Structs for binary operations f(L1,L2) to get matrix of expression - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namespace internal{ 
// L1 + L2 
struct linop_bin_add_op
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return (A.asMatrix()) + (B.asMatrix()); }
}; 
// L1 - L2 
struct linop_bin_subtract_op
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return (A.asMatrix()) - (B.asMatrix()); }
}; 
// c * L
struct scalar_left_mult_op 
{
  template<typename L2>
  auto operator()(const double& c, const L2& B) const { return  c*(B.asMatrix()); }
}; 
// -L 
struct unary_negate_op
{
  template<typename L1>
  auto operator()(const L1& B) const { return  -(B.asMatrix()); }
}; 
// composition: L1( L2( . ) )
struct linopXlinop_mult_op
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return  (A.asMatrix())*(B.asMatrix()); }
}; 

} // end namespace internal 

// Traits =====================================================================
// given a type, detect if it is derived from LinOpMixIn<> - - - - - - - - - - - - 
namespace internal{
template<typename T, typename = void>
struct is_linop_crtp_impl : std::false_type {};

template<typename T>
struct is_linop_crtp_impl<T, std::void_t<typename T::is_linop_tag>> : std::true_type {};

} // end namespace internal 

namespace traits{
template<typename T>
using is_linop_crtp = internal::is_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
} // end namespace traits 

// given a type, detect if it is a time dependent Linear Op - - - - - - - - - - - - -
namespace internal{
template<typename T, typename = void> 
struct is_time_dep_impl : public std::false_type{}; 

template<typename T> 
struct is_time_dep_impl<T,std::enable_if_t<linops::traits::is_linop_crtp<T>::value>> : public std::conditional_t<T::isTimeDep, std::true_type, std::false_type>{}; 

} // end namespace internal 

namespace traits{
template<typename T>
using is_time_dep = linops::internal::is_time_dep_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
}

// given a type, detect if it is derived from LinOpBase1D<> - - - - - - - - - - - - 
namespace internal{
template<typename T, typename = void>
struct is_1dim_linop_crtp_impl : std::false_type {};

template<typename T>
struct is_1dim_linop_crtp_impl<T, std::void_t<typename T::is_1dim_linop_tag>> : std::true_type {};

template<typename L, typename R, typename OP>
struct is_1dim_linop_crtp_impl<LinOpExpr<L,R,OP>> : std::disjunction<
  is_1dim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<L>>>,
  is_1dim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<R>>>
> {};

} // end namespace internal 

namespace traits{
template<typename T>
using is_1dim_linop_crtp = internal::is_1dim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
} // end namespace traits 

// ... from LinOpBaseXD<> - - - - - - - - - - - - 
namespace internal{
template<typename T, typename = void>
struct is_xdim_linop_crtp_impl : std::false_type {};

template<typename T>
struct is_xdim_linop_crtp_impl<T, std::void_t<typename T::is_xdim_linop_tag>> : std::true_type {};

template<typename L, typename R, typename OP>
struct is_xdim_linop_crtp_impl<LinOpExpr<L,R,OP>> : std::disjunction<
  is_xdim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<L>>>,
  is_xdim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<R>>>
> {};

} // end namespace internal 

namespace traits{
template<typename T>
using is_xdim_linop_crtp = internal::is_xdim_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
} // end namespace traits 

// Given a type, see if it is derived from CoeffMixIn - - - - - - - 
namespace internal{
template<typename T, typename = void> 
struct is_coeffop_crtp_impl : std::false_type{}; 

template<typename T>
struct is_coeffop_crtp_impl<T, std::void_t<typename T::is_coeff_flag>>: std::true_type{}; 
}

namespace traits{
template<typename T>
using is_coeffop_crtp = linops::internal::is_coeffop_crtp_impl<typename std::remove_cv_t<std::remove_reference_t<T>> >; 
} // end namespace traits

// given a type T see if it is an expression - - - - - - - - - - - - - - - - - - - - - - - - 
namespace internal{
template<typename T>
struct is_expr_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_expr_crtp_impl<LinOpExpr<L,R,OP>>: public std::true_type
{};
} // end namespace internal 

namespace traits{
template<typename T>
using is_expr_crtp = internal::is_expr_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>;
} // end namespace traits

// given a type T we may need to store it  - - - - - - - - - - - - - - - - - - - - - - - - - 
namespace traits{
// as a reference or a value in a binary expression
template<typename T, typename = void>
struct Storage 
{
  using type = T; 
}; 

template<typename LinopT>
struct Storage<LinopT, std::enable_if_t< is_linop_crtp<LinopT>::value > >
{
  using type = std::conditional_t<
    std::is_lvalue_reference<LinopT>::value && !(is_expr_crtp<LinopT>::value), // if T is an lvalue LinOpBase + not expression 
    LinopT, // store by lvalue LinopT& 
    std::remove_reference_t<LinopT> // else store rvalue by value
  >; 
}; 

} // end namespace traits 

// given type T, see if instances have .left_scalar_mult_impl(double) - - - - - - - - - - 
namespace traits{
template<typename T, typename = void>
struct supports_left_scalar_mult : public std::false_type{}; 

template<typename T>
struct supports_left_scalar_mult<T, 
  std::void_t<
    decltype(std::declval<T>().left_scalar_mult_impl( std::declval<double>() ))
  >
> : public std::true_type{}; 

} // end namespace traits 

// minimum # of doubles + return type of callable type F --------------------------------
namespace traits{
template<typename F>
class callable_traits
{
  private: 
  // ----------------------------------------------------------
  // get return type of callable type G invoked on N doubles 
  template<typename G, std::size_t N_doubles, typename... Args> 
  struct result_traits
  {
    using result_type = typename result_traits<G,N_doubles-1, double, Args...>::result_type; 
  }; 

  template<typename G, typename... Args> 
  struct result_traits<G,0,Args...>
  {
    using result_type = typename std::invoke_result<G,Args...>::type; 
  }; 

  // ---------------------------------------------------------------
  // test if callable type G can be invoked on N doubles up to max_n_args
  static constexpr std::size_t max_N_doubles = 20; 
  template<typename G, std::size_t N_doubles=max_N_doubles, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
  struct arg_traits
  {
    constexpr static bool is_callable = std::is_invocable<G,Args...>::value;

    constexpr static std::size_t num_args(){ 
      if constexpr (is_callable){
        return sizeof...(Args); 
      }
      else{
        return arg_traits<G,N_doubles-1,double, Args...>::num_args(); 
      }
    }
    
    using result_type = typename result_traits<G,num_args()>::result_type; 
  }; 

  // terminating case 
  template<typename G, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
  struct arg_traits<G,0,Args...>
  {
    constexpr static bool is_callable = std::is_invocable<G,Args...>::value; 

    constexpr static std::size_t num_args(){ 
      if constexpr (is_callable){
        return sizeof...(Args); 
      }
      else{
        static_assert(false, "maximum length of args reached"); 
      }
    } 

    using result_type = typename result_traits<F,num_args()>::result_type; 
  }; 

  // --------------------------------------------------------------------- 
  // Given a callable G, return a new type F that has constructor F( G g, double x0)
  // with an operator() that accepts # args == # args in G - 1.
  // the result is f(x1,...,xn) = g(x0,x1,...,xn) 

  template<std::size_t N_doubles, typename G, typename... Args> 
  struct BindFirst_impl : public BindFirst_impl<N_doubles-1,G,double,Args...>
  {
    using Base = BindFirst_impl<N_doubles-1,G,double,Args...>; 
    BindFirst_impl(G g, double t): Base(g,t) {}; 
    using BindFirst_impl<N_doubles-1,G,double, Args...>::operator(); 
  }; 

  template<typename G, typename... Args> 
  struct BindFirst_impl<0,G,Args...>
  {
    G func; 
    double captured; 
    BindFirst_impl(G g, double x0): func(g), captured(x0) {}; 
    double operator()(Args... args){return func(captured,args...); }; 
  }; 

  public:
  constexpr static std::size_t num_args = arg_traits<F>::num_args(); 
  using result_type = typename result_traits<F, num_args>::result_type; 
  using BindFirst_t = std::conditional_t<num_args, BindFirst_impl<num_args-1, F>, void>;
}; // end callable_traits<F> 

} // end namespace traits 

  } // end namespace linops 
} // end namespace fdm  

#endif // LinOpTraits.hpp
