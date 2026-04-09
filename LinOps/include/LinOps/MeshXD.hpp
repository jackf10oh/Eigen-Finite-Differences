// MeshXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef MESHXD_H
#define MESHXD_H 

#include<cstdint>
#include<vector>
#include<memory>
#include<algorithm>
#include<numeric>
#include<LinOps/Mesh1D.hpp>

#include<Eigen/Sparse>

namespace linops{

// forward declaration -> aliases
class MeshXD; 
using MeshXD_SPtr_t = std::shared_ptr<const MeshXD>;
using MeshXD_WPtr_t = std::weak_ptr<const MeshXD>;

class MeshXD
{
  protected:
    // Type Defs ------- 
    using StridedRef = typename Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
    using Stride_t = typename Eigen::Stride<0,Eigen::Dynamic>; 
    using StrideView_t = typename Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t>; 
    // member data --------------------------------------------------------------------------
    std::vector<Mesh1D_SPtr_t> m_mesh_vec; // dynamic array of meshes. 

  public:
    // Constructors + Destructor =============================================================== 
    // uniformly on [0,1] with 5 steps per dim, with n_dims
    MeshXD(std::size_t n_dims_init) 
      : m_mesh_vec(n_dims_init)
    {
      for(std::size_t dim_i=0; dim_i<n_dims_init; dim_i++){
        m_mesh_vec[dim_i] = linops::make_mesh(0.0, 1.0, 5); 
      }; 
    };
    
    // uniformly on [l,r] with n_steps per dim, with n_dims
    MeshXD(double left=0.0, double right=1.0, std::size_t n_steps=11,std::size_t n_dims=1)
      : m_mesh_vec(n_dims)
    {
      auto original_ptr = linops::make_mesh(left,right,n_steps); 
      for(auto& ptr : m_mesh_vec) ptr=original_ptr;  
    }
    
    // uniformly on [ L(i), R(i) ] with NSteps(i) per dim with NSteps.size() dims 
    MeshXD(const std::vector<std::pair<double,double>>& axes_vec, const std::vector<std::size_t>& nsteps_vec)
      : m_mesh_vec(axes_vec.size())
    {
      if(axes_vec.size()!=nsteps_vec.size()) throw std::invalid_argument("# of axes != # of steps given");
      auto nstep_it = nsteps_vec.begin();
      auto axes_it=axes_vec.begin();
      for(auto write=m_mesh_vec.begin(); write!=m_mesh_vec.end(); write++){
        *write = linops::make_mesh(axes_it->first, axes_it->second, *nstep_it); 
        axes_it++; 
        nstep_it++; 
      }
    }
    
    // from std::vector<> of shared_ptr<const Mesh1D>
    MeshXD(const std::vector<Mesh1D_SPtr_t>& init_vec) : m_mesh_vec(init_vec){}; 
    // from std::vector<> (RVAL) of shared_ptr<const Mesh1D>
    MeshXD(std::vector<Mesh1D_SPtr_t>&& init_vec) : m_mesh_vec(std::move(init_vec)){}; 
    
    // from shared_ptr const Mesh1D
    MeshXD(const Mesh1D_SPtr_t& m1d_init) : m_mesh_vec{m1d_init}{};  

    // Copy 
    MeshXD(const MeshXD& other)=default; 
    
    // destructors
    virtual ~MeshXD()=default;

    // Member Functions ===============================================================
    // full size of XD mesh. i.e. axis1.size() * ... * axisn.size()
    std::size_t sizes_product() const 
    {
      std::size_t p = 1; 
      for(const auto& m : m_mesh_vec) p *= m->size(); 
      return p; 
    } 

    // product of axes up to dim exclusively [first, dim)
    std::size_t sizes_middle_product(std::size_t start, std::size_t end) const 
    {
      if(start > end) throw std::invalid_argument("start index must be <= end index for middle product"); 
      if(end > m_mesh_vec.size()) throw std::invalid_argument("end index must be <= # of dims in MeshXD"); 
      std::size_t prod = 1; 
      for(std::size_t i=start; i<end; i++) prod *= m_mesh_vec[i]->size(); 
      return prod; 
    } 

    // size of a specific axis 
    std::size_t dim_size(std::size_t i) const {return m_mesh_vec.at(i)->size();} 

    // number of dimensions 
    std::size_t dims() const {return m_mesh_vec.size(); } 

    // get a specific mesh 
    Mesh1D_SPtr_t& GetMesh(std::size_t i){return m_mesh_vec[i];} 
    const Mesh1D_SPtr_t& GetMesh(std::size_t i) const {return m_mesh_vec[i];} 

    // Get a specific mesh (check index < size)
    Mesh1D_SPtr_t& GetMeshAt(std::size_t i){return m_mesh_vec.at(i);}
    const Mesh1D_SPtr_t& GetMeshAt(std::size_t i) const {return m_mesh_vec.at(i);}

    // From a VectorXd representing flattened DiscretizationXD produce list of views that "look" like 1 dimensional slices 
    std::vector<StrideView_t> OneDim_views(StridedRef vec, std::size_t ith_dim=0) const 
    {
      // # of entries in vec must be == to product of mesh1D sizes
      if(vec.size() != sizes_product()) throw std::runtime_error("DiscretizationXD # of entries must be == to product of sizes in MeshXD"); 
      // i has to be one of the dimensions of DiscretizationXD 
      if(ith_dim >= dims()) throw std::runtime_error("MeshXD::OneDim_views(i) i must be < MeshXD.dims().");

      std::size_t ith_dim_size = dim_size(ith_dim); 
      std::size_t num_copies = sizes_product() / ith_dim_size; 
      std::size_t mod = sizes_middle_product(0, ith_dim); 
      std::size_t scale = mod * ith_dim_size; 

      std::vector<StrideView_t> result; 
      result.reserve(dim_size(ith_dim)); 

      Stride_t stride(0,mod); 

      // iterate through the copies 
      for(std::size_t n=0; n<num_copies; n++)
      {
        // offset from start of current copy 
        std::size_t offset = (mod ? n % mod : n) + (scale * (n/mod));  
        // begin data ptr of copy  
        auto begin = vec.data()+offset;  

        // MemView of current copy 
        result.emplace_back(begin, ith_dim_size, stride);
      }
      return result; 
    } 
};

template<typename... Args> 
auto make_meshXD(Args... args)
{
  return std::make_shared<const MeshXD>(args...); 
}

// just return the input args if called on a shared_ptr<MeshXD>
auto make_meshXD(const std::shared_ptr<const MeshXD>& other)
{
  return other; 
}

} // end namespace linops 

#endif // MeshXD.hpp