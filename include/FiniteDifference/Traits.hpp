// Traits.hpp
//
// traits used by all subdirectory of fdm library 
//
// JAF 4/14/2026 

#ifndef FDM_TRAITS_H
#define FDM_TRAITS_H

#include<cstdint>
#include<string>
#include<type_traits> // decay_t 
#include<complex>
#include<Eigen/Core>
#include<Eigen/src/Core/util/Macros.h>
#include<Eigen/src/Core/util/Constants.H> 
#include<Eigen/src/Core/util/ForwardDeclarations.h>  // CwiseUnaryOp, CwiseBinaryOp, 
#include<Eigen/src/Core/EigenBase.h> 
#include<Eigen/src/SparseCore/SparseUtil.h> // forward declares SparseMatrix<...> 
#include<Eigen/src/SparseCore/CompressedStorage.h>
#include<Eigen/src/SparseCore/SparseCompressedBase.h>
// #include<Eigen/SparseCore> can't include before plugin macro takes effect! 
#include "Types.hpp"

namespace fdm{ 
namespace internal{

// traits around a callable F ----------------------------------------------
template<typename F>
class callable_traits
{
  private: 
  // get return type of callable type G invoked on N scalars 
  template<typename G, std::size_t numScalars, typename... Args> 
  struct result_traits
  {
    using result_type = typename result_traits<G,numScalars-1, fdm::Scalar, Args...>::result_type; 
  }; 

  template<typename G, typename... Args> 
  struct result_traits<G,0,Args...>
  {
    using result_type = typename std::invoke_result<G,Args...>::type; 
  }; 

  // ---------------------------------------------------------------
  // test if callable type G can be invoked on N scalars up to max_n_args
  static constexpr std::size_t numScalarsMax = 20; 
  template<typename G, std::size_t numScalars=numScalarsMax, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
  struct arg_traits
  {
    constexpr static bool is_callable = std::is_invocable<G,Args...>::value;

    constexpr static std::size_t num_args(){ 
      if constexpr (is_callable){
        return sizeof...(Args); 
      }
      else{
        return arg_traits<G,numScalars-1,fdm::Scalar, Args...>::num_args(); 
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
  // Given a callable G, return a new type F that has constructor F( G g, Scalar x0)
  // with an operator() that accepts # args == # args in G - 1.
  // the result is f(x1,...,xn) = g(x0,x1,...,xn) 

  template<std::size_t numScalars, typename G, typename... Args> 
  struct BindFirst_impl : public BindFirst_impl<numScalars-1,G,fdm::Scalar,Args...>
  {
    using Base = BindFirst_impl<numScalars-1,G,fdm::Scalar,Args...>; 
    BindFirst_impl(G g, fdm::Scalar t): Base(g,t) {}; 
    using BindFirst_impl<numScalars-1,G,fdm::Scalar, Args...>::operator(); 
  }; 

  template<typename G, typename... Args> 
  struct BindFirst_impl<0,G,Args...>
  {
    const G func; 
    fdm::Scalar captured; 
    BindFirst_impl(G g, fdm::Scalar x0): func(g), captured(x0) {}; 
    fdm::Scalar operator()(Args... args) const {return func(captured,args...); }; 
  }; 

  public:
  constexpr static std::size_t num_args = arg_traits<F>::num_args(); 
  using result_type = typename result_traits<F, num_args>::result_type; 
  using BindFirst = std::conditional_t<num_args, BindFirst_impl<num_args-1, F>, void>;
}; // end callable_traits

} // end namespace internal 
} // end namespace fdm

#endif // Traits.hpp 