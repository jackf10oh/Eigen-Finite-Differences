// VectorXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef VECTORXD_H
#define VECTORXD_H

#include<vector>
#include<Eigen/Core>
#include "MeshXD.hpp"
#include<LinOps/LinOpTraits.hpp> // callable_traits<T> 

namespace linops{

class VectorXD
{
  private:
    // private typedefs ==================================================== 
    typedef Eigen::Stride<0,Eigen::Dynamic> Stride_t; 
    typedef Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t> StrideView_t; 

    // member data ==========================================================
    Eigen::VectorXd m_vals; // flattened array of values 
    MeshXD_WPtr_t m_mesh_ptr; 

  public:
    // Constructors / Destructors  ==========================================================
    // Default ----------------------------------------------- 
    VectorXD()=default; 
    // from size of m_vals. assume just 1 dimension ----------
    VectorXD(std::size_t size_init): m_mesh_ptr(), m_vals(size_init){}; 
    // from a MeshXDPtr --------------------------------------
    VectorXD(const MeshXD_SPtr_t& mesh_init) 
      : m_mesh_ptr(mesh_init), m_vals(mesh_init->sizes_product())
    {}
    // copy --------------------------------------
    VectorXD(const VectorXD& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {}; 

    // copy from Eigen::VectorXd
    VectorXD(const Eigen::VectorXd& other, MeshXD_WPtr_t mesh_init = MeshXD_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(other)
    {}; 

    // move ---------------------------------------- 
    VectorXD(VectorXD&& other)
    : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {}; 
    
    // move from Eigen::VectorXD
    VectorXD(Eigen::VectorXd&& other, MeshXD_WPtr_t mesh_init = MeshXD_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(std::move(other))
    {}; 
    // destructor 
    ~VectorXD()=default; 

    // member functions ==========================================================
    // get underlying values 
    Eigen::VectorXd& values(){return m_vals; }
    const Eigen::VectorXd& values() const {return m_vals; } 

    // Give a list of Eigen::Map<>. each Map looks like a Vector1D on a Mesh1d  
    std::vector<StrideView_t> OneDim_views(std::size_t ith_dim=0)
    {
      return m_mesh_ptr.lock()->OneDim_views(m_vals, ith_dim); 
    }

    // get underlying MeshXD_SPtr_t 
    MeshXD_SPtr_t get_meshxd() const {return m_mesh_ptr.lock();}; 

    // number of dimensions 
    std::size_t dims() const{ return m_mesh_ptr.lock()->dims(); };
    // size of ith dimension  
    std::size_t dim_size(std::size_t i) const {return m_mesh_ptr.lock()->dim_size(i); }
    
    // product of all dimensions' sizes 
    std::size_t sizes_product() const { return m_vals.size(); } 
    // product of dimensions in [start,end)
    std::size_t sizes_middle_product(std::size_t start, std::size_t end){
      return m_mesh_ptr.lock()->sizes_middle_product(start,end); 
    }

    // store a new MeshXD_WPtr_t
    VectorXD& set_mesh(MeshXD_WPtr_t m){ m_mesh_ptr = m; return *this; } 
    // set discretization to same size as meshxd's sizes_product
    VectorXD& resize(const MeshXD_SPtr_t& m) { 
      m_mesh_ptr=m; 
      m_vals.conservativeResize(m->sizes_product()); 
      return *this; 
    }
       
    // Operators ----------------------------------------------------
    VectorXD& operator=(const VectorXD& other) = default;
    VectorXD& operator=(VectorXD&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      // other.m_mesh_ptr = nullptr; 
      m_vals = std::move(other.m_vals); 
      return *this;
    }; 
}; 

// set vector to match a mesh size and set it constant 
linops::VectorXD make_Discretization(const MeshXD_SPtr_t& m, double val){ 
  linops::VectorXD result(m);
  result.values().setConstant(val); 
  return result;
} 

// set discretizations values according to callable type F
template<
typename F,
typename = std::enable_if_t<
  !std::is_arithmetic_v<std::remove_reference_t<std::remove_cv_t<F>>>
  >
>
linops::VectorXD make_Discretization(const MeshXD_SPtr_t& m, F func)
{
  // assert func returns double 
  static_assert(std::is_same<typename traits::callable_traits<F>::result_type, double>::value, "static assert error: callable type F must return a double"); 
  
  // check there are enough dimensions to use callable type F
  static constexpr std::size_t num_args = traits::callable_traits<F>::num_args; 
  if(m->dims() < num_args) 
    throw std::invalid_argument(
      "# dims of MeshXD_SPtr_t must be >= # args in callable F"); 
  
  // result returns by make_Discretization() 
  linops::VectorXD result(m); 

  if constexpr(num_args == 0){
    result.values().setConstant( func() );
    return result; 
  }

  // stores sizes_middle_product(0,i) for i=0,...,num_args
  std::array<std::size_t, num_args> cumulative_prod_arr; 
  for(std::size_t d=0; d<num_args; d++) cumulative_prod_arr[d] = m->sizes_middle_product(0, d); 
  
  // iterate through flattened first layer 
  std::array<double, num_args> args_arr; // stores n args for std::apply(func,args_arr) later
  std::size_t flat_end = m->sizes_middle_product(0,num_args);   
  for(std::size_t flat_i=0; flat_i<flat_end; flat_i++){
    // fill args_arr 
    for(std::size_t d=0; d < num_args; d++){
      std::size_t stride = cumulative_prod_arr[d]; 
      std::size_t index = (flat_i / stride) % m->dim_size(d); 
      args_arr[d] = (*m->GetMesh(d))[index]; 
    }
    // set all values according to func(x0, x1, ... , xn) 
    result.values()[flat_i] = std::apply(func, args_arr);  
  } // end first layer 

  // if we need to copy into more layers 
  if(flat_end != result.values().size()){
    // ither through views
    std::size_t n_layers = m->sizes_product() / flat_end; // # of times to copy/paste first [0,flat_end) vales
    for(std::size_t ith_view=0; ith_view<flat_end; ith_view++){
      // fill in all layers with first layers value 
      for(std::size_t layer=1; layer<n_layers; layer++){
        result.values()[ith_view + flat_end*layer] = result.values()[ith_view]; 
      } // end for loop through values of ith layer 
    } // end for loop through layers  
  } // end if 

  // all values filled 
  return result;
}

} // end namespace linops 

#endif // VectorXD.hpp