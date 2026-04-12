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

#include "Mesh1D.hpp"

namespace fdm{

// forward declaration -> aliases
class MeshXD; 
using SharedConstMeshXD = std::shared_ptr<const MeshXD>;
using WeakConstMeshXD = std::weak_ptr<const MeshXD>;

class MeshXD
{
  private:
    // Types -----------------
    using StridedRef = typename Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>; 
    using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
    using StrideView =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride>; 

    // member data --------------------------------------------------------------------------
    std::vector<SharedConstMesh1D> m_mesh_vec; // dynamic array of meshes. 

  public:
    // Constructors + Destructor =============================================================== 
    
    // uniformly on [0,1] with 5 steps per dim, with n_dims
    MeshXD(std::size_t n_dims_init) 
      : m_mesh_vec(n_dims_init)
    {
      for(std::size_t dim_i=0; dim_i<n_dims_init; dim_i++){
        m_mesh_vec[dim_i] = fdm::make_Mesh1D(0.0, 1.0, 5); 
      }; 
    };
    
    // uniformly on [l,r] with n_steps per dim, with n_dims
    MeshXD(double left=0.0, double right=1.0, std::size_t n_steps=11,std::size_t n_dims=1)
      : m_mesh_vec(n_dims)
    {
      auto original_ptr = fdm::make_Mesh1D(left,right,n_steps); 
      for(auto& ptr : m_mesh_vec) ptr=original_ptr;  
    }
    
    // uniformly on [ L(i), R(i) ] with NSteps(i) per dim with NSteps.size() numDims 
    MeshXD(const std::vector<std::pair<double,double>>& axes_vec, const std::vector<std::size_t>& nsteps_vec)
      : m_mesh_vec(axes_vec.size())
    {
      if(axes_vec.size()!=nsteps_vec.size()) throw std::invalid_argument("# of axes != # of steps given");
      auto nstep_it = nsteps_vec.begin();
      auto axes_it=axes_vec.begin();
      for(auto write=m_mesh_vec.begin(); write!=m_mesh_vec.end(); write++){
        *write = fdm::make_Mesh1D(axes_it->first, axes_it->second, *nstep_it); 
        axes_it++; 
        nstep_it++; 
      }
    }
    
    // from std::vector<> of shared_ptr<const Mesh1D>
    MeshXD(const std::vector<SharedConstMesh1D>& init_vec) : m_mesh_vec(init_vec){}; 
    
    // from std::vector<> (RVAL) of shared_ptr<const Mesh1D>
    MeshXD(std::vector<SharedConstMesh1D>&& init_vec) : m_mesh_vec(std::move(init_vec)){}; 
    
    // from shared_ptr const Mesh1D
    MeshXD(const SharedConstMesh1D& m1d_init) : m_mesh_vec{m1d_init}{};  

    // Copy 
    MeshXD(const MeshXD& other)=default; 
    
    // destructors
    ~MeshXD()=default;

    // Member Functions ===============================================================
    // full size of XD mesh. i.e. axis1.size() * ... * axisn.size()
    std::size_t sizesProduct() const 
    {
      std::size_t p = 1; 
      for(const auto& m : m_mesh_vec) p *= m->size(); 
      return p; 
    } 

    // product of axes up to dim exclusively [first, dim)
    std::size_t sizesMiddleProduct(std::size_t start, std::size_t end) const 
    {
      if(start > end) throw std::invalid_argument("start index must be <= end index for middle product"); 
      if(end > m_mesh_vec.size()) throw std::invalid_argument("end index must be <= # of numDims in MeshXD"); 
      std::size_t prod = 1; 
      for(std::size_t i=start; i<end; i++) prod *= m_mesh_vec[i]->size(); 
      return prod; 
    } 

    // size of a specific axis 
    std::size_t sizeOfDim(std::size_t i) const {return m_mesh_vec.at(i)->size();} 

    // number of dimensions 
    std::size_t numDims() const {return m_mesh_vec.size(); } 

    // get a specific mesh 
    SharedConstMesh1D& getMesh1D(std::size_t i){return m_mesh_vec[i];} 
    const SharedConstMesh1D& getMesh1D(std::size_t i) const {return m_mesh_vec[i];} 

    // Get a specific mesh (check index < size)
    SharedConstMesh1D& getMesh1DSafe(std::size_t i){return m_mesh_vec.at(i);}
    const SharedConstMesh1D& getMesh1DSafe(std::size_t i) const {return m_mesh_vec.at(i);}

    // From a VectorXd representing flattened DiscretizationXD produce list of views that "look" like 1 dimensional slices 
    std::vector<StrideView> makeOneDimViews(StridedRef vec, std::size_t ith_dim=0) const 
    {
      // # of entries in vec must be == to product of mesh1D sizes
      if(vec.size() != sizesProduct()) throw std::runtime_error("DiscretizationXD # of entries must be == to product of sizes in MeshXD"); 
      // i has to be one of the dimensions of DiscretizationXD 
      if(ith_dim >= numDims()) throw std::runtime_error("MeshXD::OneDim_views(i) i must be < MeshXD.numDims().");

      std::size_t ith_dim_size = sizeOfDim(ith_dim); 
      std::size_t num_copies = sizesProduct() / ith_dim_size; 
      std::size_t mod = sizesMiddleProduct(0, ith_dim); 
      std::size_t scale = mod * ith_dim_size; 

      std::vector<StrideView> result; 
      result.reserve(sizeOfDim(ith_dim)); 

      Stride stride(0,mod); 

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
auto make_MeshXD(Args... args)
{
  return std::make_shared<const MeshXD>(args...); 
}

// just return the input args if called on a shared_ptr<MeshXD>
const auto& make_MeshXD(const std::shared_ptr<const MeshXD>& other)
{
  return other; 
}

} // end namespace fdm 

#endif // MeshXD.hpp