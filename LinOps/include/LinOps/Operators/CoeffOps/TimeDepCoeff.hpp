// TimeDepCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include<Eigen/Core> // VectorXd
#include<Utilities/SparseDiagExpr.hpp> // SparseDiag

#include "../../CoeffOpMixIn.hpp"
#include "../../LinOpExpr.hpp" 
#include "../../LinOpTraits.hpp" // callable_traits<F>
#include "../../Vector.hpp"
#include "../../VectorXD.hpp"

namespace linops{

namespace internal{

template<std::size_t N, typename = std::enable_if_t<(N>0), void> >
struct TimeDepMemberData
{
  Eigen::VectorXd m_diag_vals = Eigen::VectorXd(0);  
  std::size_t m_prod_after = 1; 
};

template<>
struct TimeDepMemberData< 1, void>
{
  double m_scalar_val; 
};

} // end namespace internal 

template<typename FUNC_STORAGE_T>
class TimeDepCoeff :
  public CoeffOpMixIn<TimeDepCoeff<FUNC_STORAGE_T>>, 
  public LinOpBase1D<TimeDepCoeff<FUNC_STORAGE_T>>, 
  public LinOpBaseXD<TimeDepCoeff<FUNC_STORAGE_T>>, 
  private internal::TimeDepMemberData<traits::callable_traits<FUNC_STORAGE_T>::num_args>
{
  public:
    // flags as time dependent. 
    constexpr static bool is_time_dep_flag = true; 

  private:
    constexpr static std::size_t N = traits::callable_traits<FUNC_STORAGE_T>::num_args; 
    // Member Data -----------------------------------
    // the callable object 
    FUNC_STORAGE_T m_function;  
    // TimeDepCoeff stores an owned copy of std::shared_ptr<MeshXD> for first N dimensions so SetTime_impl(t) doesn't incur reference counting.  
    MeshXD_SPtr_t m_owned_subdim_mesh_ptr; 
    // point to full XD mesh. 
    std::variant<Mesh1D_WPtr_t, MeshXD_WPtr_t> m_mesh_ptr; 

  public:
    // Constructors + Destructor ==========================================================
    TimeDepCoeff()=delete; // no default constructor

    // from callable + mesh1d
    TimeDepCoeff(FUNC_STORAGE_T f_init, const linops::Mesh1D_SPtr_t& m = nullptr)
      : m_function(f_init)
    {
      constexpr bool returns_double = std::is_same<double, typename linops::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value; 
      constexpr bool args_more_than_zero = (traits::callable_traits<FUNC_STORAGE_T>::num_args > 0);
      static_assert(args_more_than_zero && returns_double,"Error constructing coeff: F doesn't return double or have args > 0");  
      if(m)
      {
        set_mesh(m);
        // SetTime_impl( this->m_current_time ); // move into separate constructor...
      }
    }
    // from callable + meshxd 
    TimeDepCoeff(FUNC_STORAGE_T f_init, const linops::MeshXD_SPtr_t& m)
      : m_function(f_init)
    {
      constexpr bool returns_double = std::is_same<double, typename linops::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value; 
      constexpr bool args_more_than_zero = (traits::callable_traits<FUNC_STORAGE_T>::num_args > 0);
      static_assert(args_more_than_zero && returns_double,"Error constructing coeff: F doesn't return double or have args > 0");  
      if(m)
      {
        set_mesh(m);
        // SetTime_impl( this->m_current_time ); // move into separate constructor...
      }
    }
    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other)=default; 
    
    // destructor 
    ~TimeDepCoeff()=default;

    // Member Funcs =============================================================
    // Matrix getters 
    auto asMatrix()
    {
      if constexpr(N > 1){
        return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(this->m_diag_vals, this->m_prod_after); 
      }
      else if constexpr(N == 1){
        return this->m_scalar_val; 
      }
      else{
        static_assert(false, "must have num_args >= 1 in TimeDep"); 
      }
    };
    auto asMatrix() const
    {
      if constexpr(N > 1){
        return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(this->m_diag_vals, this->m_prod_after); 
      }
      else if constexpr(N == 1){
        return this->m_scalar_val; 
      }
      else{
        static_assert(false, "must have num_args >= 1 in TimeDep"); 
      }
    };

    // return weak_ptr of Mesh1D pointed to
    Mesh1D_WPtr_t get_weak_mesh1d() const { 
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr); 
      return Mesh1D_WPtr_t{}; 
    }

    // return weak_ptr of MeshXD pointed to
    MeshXD_WPtr_t get_weak_meshxd() const { 
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr); 
      return MeshXD_WPtr_t{}; 
    }

    // return Mesh1D pointed to 
    Mesh1D_SPtr_t get_mesh1d() const {
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr).lock(); 
      return Mesh1D_SPtr_t{};
    } 

    // return MeshXD pointed to 
    MeshXD_SPtr_t get_meshxd() const {
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr).lock(); 
      return MeshXD_SPtr_t{};
    } 

    // updated m_mesh_ptr (1D)
    void set_mesh(const Mesh1D_SPtr_t& m){
      this->m_mesh_ptr = m; 
      m_owned_subdim_mesh_ptr = make_meshXD(m); 
      if constexpr(N > 1 ){
        this->m_prod_after = 1; 
      }
    }

    // updated m_mesh_ptr (XD)
    void set_mesh(const MeshXD_SPtr_t& m){
      this->m_mesh_ptr = m; 

      // Take ownership of first N-1 mesh1d. throws error if N-1 > # of dims in m
      std::vector<Mesh1D_SPtr_t> v; 
      for(auto i=0; i< N-1; i++) v.emplace_back(m->GetMeshAt(i));  
      m_owned_subdim_mesh_ptr = make_meshXD(std::move(v)); 
      if constexpr(N > 1){
        this->m_prod_after = m->sizes_product() / m_owned_subdim_mesh_ptr->sizes_product();  
      }
    }

    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){

      if constexpr(N == 1){
        this->m_scalar_val = m_function(t); 
      }
      else if constexpr(N > 1){
        using Bind_t = typename traits::callable_traits<FUNC_STORAGE_T>::BindFirst_t; 
        Bind_t binded(m_function, t); 
        this->m_diag_vals = linops::make_Discretization(m_owned_subdim_mesh_ptr,binded).values(); 
      }
      else{
        static_assert(false, "must use # of args >= 1 in Func for TimeDepCoeff"); 
      }
    };
}; 

} // end namespace linops 

#endif // TimeDepCoeff.hpp 