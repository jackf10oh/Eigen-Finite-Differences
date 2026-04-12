// LinOpMixIn.hpp 
//
// CRTP mixin class for both LinOpBase1D + XD
//
// JAF 1/29/2026

#ifndef LINOPMIXIN_H
#define LINOPMIXIN_H

#include "LinOpTraits.hpp"
#include "LinOpExpr.hpp"

namespace fdm{
  namespace linops{

namespace internal{
// Leafs of an expression hold double m_current_time
template<typename T, typename = void>
struct LinOpMixInData{
  double m_current_time=0.0; 
};

// expression itself holds no data 
template<typename T>
struct LinOpMixInData<T, std::enable_if_t< linops::traits::is_expr_crtp<T>::value >>
{};
} // end namespace internal 

template<typename Derived>
class LinOpMixIn : private internal::LinOpMixInData<Derived>
{
  // optional macro to add features into both LinOpBase1D<> + XD<>
  #ifndef LINOP_PLUGIN
  #else
  #include LINOP_PLUGIN
  #endif

  public:
    // Type Defs -------------------------------------
    // to tell if a class derived from LinOpMixIn<> 
    typedef struct{} is_linop_tag; 
    // so LinOpMixIn can access most derived class 
    using DerivedT = Derived;

    // Member Data -------------------
    // detects if Derived implements its own setTime_impl. 
    static constexpr bool isTimeDep = (&DerivedT::setTime_impl != &LinOpMixIn::setTime_impl); 
    // double m_current_time inherited in LinOpMixInData<>

  public:

    // Member Funcs ==========================================================

    // Get turrent time ----------------------------------
    const double& getTime() const 
    {
      // if we are an expression 
      if constexpr(traits::is_expr_crtp<typename Derived::DerivedT>::value){
        const auto& expr = static_cast<const typename Derived::DerivedT&>(*this);
        if constexpr(linops::traits::is_linop_crtp<typename Derived::DerivedT::LStorage>::value){
          return expr.getLhs().getTime();
        }
        else if constexpr(linops::traits::is_linop_crtp<typename Derived::DerivedT::RStorage>::value){
          return expr.Rhs().getTime();
        }
      }
      // otherwise 
      else{
        return this->m_current_time; 
      }
    } 

    // Set the current time of the operator. ------------------
    void setTime(double t)
    {
      // if we are an expression 
      if constexpr (traits::is_expr_crtp<typename Derived::DerivedT>::value)
      {
        auto& expr = static_cast<typename Derived::DerivedT&>(*this);
        // if LHS of expr is LinOp
        if constexpr(linops::traits::is_linop_crtp<typename Derived::DerivedT::LStorage>::value) 
        {
          // LHS sets time 
          expr.getLhs().setTime(t);
        }
        // if RHS of expr is LinOp
        if constexpr(linops::traits::is_linop_crtp<typename Derived::DerivedT::RStorage>::value)
        {
          // RHS sets time 
          expr.getRhs().setTime(t);
        } 
      }
      // we aren't an expression call the actual implementor 
      else  
      {
        // store new time.   
        this->m_current_time = t;
        // static_cast<typename Derived::DerivedT*>(this)->setTime_impl(t);
        Accessor::call_setTime_impl(*this, t); 
      } 
    }
    void setTime_impl(double t){ /* do nothing by default ...*/}; 

    // convert a linear operator to its matrix form -----------------------
    decltype(auto) asMatrix() const
    {
      return static_cast<const typename Derived::DerivedT*>(this)->asMatrix(); 
    };

    // Composition of Linear Ops L1(L2( . )) (lval) ---------------------------------------------------
    template<typename R> 
    auto compose(R&& InnerOp) &
    {
      static_assert(traits::is_linop_crtp<R>::value, "compose only works on other linear operators!"); 
      using Lhs = DerivedT&; 
      using Rhs = std::remove_reference_t<R>;
      return LinOpExpr<Lhs, Rhs, internal::linopXlinop_mult_op>(
      std::forward<Lhs>(static_cast<Lhs>(*this)), // lhs
      std::forward<R>(InnerOp) // rhs 
      );
    }; // end .compose(other) & lvalue overload 

    // // composition of linear of L1(L2( . )) (rval)
    template<typename R>
    auto compose(R&& InnerOp) && 
    {
      static_assert(traits::is_linop_crtp<R>::value, "compose only works on other linear operators!"); 
      using Lhs = DerivedT&&; 
      using Rhs = std::remove_reference_t<R>;
      return LinOpExpr<Lhs, Rhs, internal::linopXlinop_mult_op>(
      std::forward<Lhs>(static_cast<Lhs>(*this)), // lhs
      std::forward<R>(InnerOp) // rhs 
      ); 
    }; // end .compose(other) && rvalue overload  

  private:
    // not accessibles ==============================================================================    
    // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    template<typename Scalar>
    auto left_scalar_mult_impl(Scalar&& c) & {
      return LinOpExpr<Scalar, typename Derived::DerivedT&, internal::scalar_left_mult_op>(
        std::forward<Scalar>(c), // lhs scalar
        static_cast<typename Derived::DerivedT&>(*this) // rhs 
      );
    }
    
    // Left multiply by a scalar: i.e. c*L (rval)
    template<typename Scalar>
    auto left_scalar_mult_impl(Scalar&& c) && {
      return LinOpExpr<Scalar, typename Derived::DerivedT&&, internal::scalar_left_mult_op>(
        std::forward<Scalar>(c), // lhs scalar
        static_cast<typename Derived::DerivedT&&>(*this) // rhs 
      );
    }

  public:
    // Operators ================================================================================ 
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename R, typename = std::enable_if_t<traits::is_linop_crtp<R>::value>>
    auto operator+(R&& rhs) & {
        return LinOpExpr<typename Derived::DerivedT&, R, internal::linop_bin_add_op>(
        static_cast<typename Derived::DerivedT&>(*this),  // lhs
        std::forward<R>(rhs) // rhs 
      );
    }
    
    // L1 + L2 (rval)  
    template<typename R, typename = std::enable_if_t<traits::is_linop_crtp<R>::value>>
    auto operator+(R&& rhs) && {
        return LinOpExpr<typename Derived::DerivedT&&, R, internal::linop_bin_add_op>(
        static_cast<typename Derived::DerivedT&&>(*this), // lhs 
        std::forward<R>(rhs) // rhs
      );
    }
    
    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename R, typename = std::enable_if_t<traits::is_linop_crtp<R>::value>>
    auto operator-(R&& rhs) & {
        return LinOpExpr<typename Derived::DerivedT&, R, internal::linop_bin_subtract_op>(
        static_cast<typename Derived::DerivedT&>(*this), //lhs
        std::forward<R>(rhs) //rhs
      );
    }
    
    // L1 - L2 (rval)  
    template<typename R, typename = std::enable_if_t<traits::is_linop_crtp<R>::value>>
    auto operator-(R&& rhs) && {
        return LinOpExpr<typename Derived::DerivedT&&, R, internal::linop_bin_subtract_op>(
        static_cast<typename Derived::DerivedT&&>(*this), // lhs
        std::forward<R>(rhs) // rhs 
      );
    }
    
    // friend declare c * L (lval + rval) -------------------------------------------
    template<typename Scalar, typename R, typename>
    friend auto operator*(Scalar&& scalar, R&& rhs); 
    template<typename T>
    friend struct traits::supports_left_scalar_mult;

    // unary operator-() (lval) ---------------------------------------------- 
    auto operator-() & {
        return LinOpExpr<typename Derived::DerivedT&, void, internal::unary_negate_op>(
        static_cast<typename Derived::DerivedT&>(*this) // lhs 
      );
    }
    // unary operator-() (rval) 
    auto operator-() && {
        return LinOpExpr<typename Derived::DerivedT&&, void, internal::unary_negate_op>(
        static_cast<typename Derived::DerivedT&&>(*this) // lhs
      );
    }

  private:
    struct Accessor : public DerivedT
    {
      static void call_setTime_impl(LinOpMixIn& self, double t){ static_cast<Accessor&>(self).setTime_impl(t); }
    };

}; // end LinOpMixIn<> 

// operator*(c,L) outside of class ....
template<
  typename Scalar, 
  typename R, 
  typename = std::enable_if_t<
    std::conjunction_v<
    traits::is_linop_crtp<R>, 
    traits::supports_left_scalar_mult<R>, 
    std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<Scalar>>>
    > // end conjuntion 
  > // end enable_if
>
auto operator*(Scalar&& c, R&& rhs){
  return std::forward<R>(rhs).left_scalar_mult_impl( std::forward<Scalar>(c) ); 
}

  } // end namespace linops 
} // end namespace fdm  

#endif // LinOpMixIn.hpp 
