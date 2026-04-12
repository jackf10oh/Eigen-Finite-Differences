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
#include "LinOps/LinOpTraits.hpp" // callable_traits<T> 

namespace fdm{

class VectorXD
{
  private:
    // member data ---------------------------------------------------------------
    Eigen::VectorXd m_vals; // flattened array of values 
    WeakConstMeshXD m_mesh_ptr; 

  public:
    // Constructors / Destructors  ==========================================================
    
    // default 
    VectorXD()=default; 
    
    // from size of m_vals. assume just 1 dimension
    VectorXD(std::size_t size_init): m_mesh_ptr(), m_vals(size_init){}; 
    
    // from a MeshXDPtr 
    VectorXD(const SharedConstMeshXD& mesh_init) 
      : m_mesh_ptr(mesh_init), m_vals(mesh_init->sizesProduct())
    {}
    
    // copy 
    VectorXD(const VectorXD& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {}; 

    // copy from Eigen::VectorXd
    VectorXD(const Eigen::VectorXd& other, WeakConstMeshXD mesh_init = WeakConstMeshXD{})
      : m_mesh_ptr(mesh_init), m_vals(other)
    {}; 

    // move  
    VectorXD(VectorXD&& other)
    : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {}; 
    
    // move from Eigen::VectorXD
    VectorXD(Eigen::VectorXd&& other, WeakConstMeshXD mesh_init = WeakConstMeshXD{})
      : m_mesh_ptr(mesh_init), m_vals(std::move(other))
    {}; 
    // destructor 
    ~VectorXD()=default; 

    // member functions ==========================================================
    
    // get underlying values 
    Eigen::VectorXd& values(){return m_vals; }
    const Eigen::VectorXd& values() const {return m_vals; } 

    // Give a list of Eigen::Map<>. each Map looks like a Vector1D on a Mesh1d  
    auto getOneDimViews(std::size_t ith_dim=0)
    {
      return m_mesh_ptr.lock()->makeOneDimViews(m_vals, ith_dim); 
    }

    // get underlying SharedConstMeshXD 
    SharedConstMeshXD get_meshxd() const {return m_mesh_ptr.lock();}; 

    // number of dimensions 
    std::size_t numDims() const{ return m_mesh_ptr.lock()->numDims(); };
    
    // size of ith dimension  
    std::size_t sizeOfDim(std::size_t i) const {return m_mesh_ptr.lock()->sizeOfDim(i); }
    
    // product of all dimensions' sizes 
    std::size_t sizesProduct() const { return m_vals.size(); } 
    
    // product of dimensions in [start,end)
    std::size_t sizesMiddleProduct(std::size_t start, std::size_t end){
      return m_mesh_ptr.lock()->sizesMiddleProduct(start,end); 
    }

    // store a new WeakConstMeshXD
    VectorXD& set_mesh(WeakConstMeshXD m){ m_mesh_ptr = m; return *this; } 
    // set discretization to same size as meshxd's sizesProduct
    VectorXD& resize(const SharedConstMeshXD& m) { 
      m_mesh_ptr=m; 
      m_vals.conservativeResize(m->sizesProduct()); 
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
fdm::VectorXD make_Discretization(const SharedConstMeshXD& m, double val){ 
  fdm::VectorXD result(m);
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
fdm::VectorXD make_Discretization(const SharedConstMeshXD& m, F func)
{
  // assert func returns double 
  static_assert(std::is_same<typename linops::traits::callable_traits<F>::result_type, double>::value, "static assert error: callable type F must return a double"); 
  
  // check there are enough dimensions to use callable type F
  static constexpr std::size_t num_args = linops::traits::callable_traits<F>::num_args; 
  if(m->numDims() < num_args) 
    throw std::invalid_argument(
      "# numDims of SharedConstMeshXD must be >= # args in callable F"); 
  
  // result returns by make_Discretization() 
  fdm::VectorXD result(m); 

  if constexpr(num_args == 0){
    result.values().setConstant( func() );
    return result; 
  }

  // stores sizesMiddleProduct(0,i) for i=0,...,num_args
  std::array<std::size_t, num_args> cumulative_prod_arr; 
  for(std::size_t d=0; d<num_args; d++) cumulative_prod_arr[d] = m->sizesMiddleProduct(0, d); 
  
  // iterate through flattened first layer 
  std::array<double, num_args> args_arr; // stores n args for std::apply(func,args_arr) later
  std::size_t flat_end = m->sizesMiddleProduct(0,num_args);   
  for(std::size_t flat_i=0; flat_i<flat_end; flat_i++){
    // fill args_arr 
    for(std::size_t d=0; d < num_args; d++){
      std::size_t stride = cumulative_prod_arr[d]; 
      std::size_t index = (flat_i / stride) % m->sizeOfDim(d); 
      args_arr[d] = (*m->getMesh1D(d))[index]; 
    }
    // set all values according to func(x0, x1, ... , xn) 
    result.values()[flat_i] = std::apply(func, args_arr);  
  } // end first layer 

  // if we need to copy into more layers 
  if(flat_end != result.values().size()){
    // ither through views
    std::size_t n_layers = m->sizesProduct() / flat_end; // # of times to copy/paste first [0,flat_end) vales
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

} // end namespace fdm 

#endif // VectorXD.hpp