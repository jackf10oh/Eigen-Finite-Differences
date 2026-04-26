// Mesh.hpp 
//
// list of Eigen::VectorXd representing
// a disretization of space in any # of dimensions 
// 
// JAF 4/13/2026 

#ifndef ANYDIMENSIONALMESH_H
#define ANYDIMENSIONALMESH_H 

#include<vector> // std::vector
#include<memory> // std::shared_ptr
#include<Eigen/Core> // Eigen::VectorXd 
#include<iostream> // TODO remove this 
// #include "LinOps/LinOpTraits.hpp"

namespace fdm {

// forward declaration + alias
class Mesh; 
using SharedMesh = std::shared_ptr<Mesh>; 
using SharedConstMesh = std::shared_ptr<const Mesh>; 

class Mesh : public std::enable_shared_from_this<Mesh>
{
  private:
    // Type Defs ----------------------- 
    using Stride =  Eigen::Stride<0,Eigen::Dynamic>; 
    using StrideView =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride>;
    // member data ----------------------------
    static constexpr std::size_t numDimsMax = 5; // fixed maximum.
    typename std::array<Eigen::VectorXd, numDimsMax> m_mesh_arr; // fixed size array of 1d axes.
    std::size_t m_size; // runtime size that counts how many dims are used. 

  public:
    // Constructors + Destructor =============================================================== 

    // from number of axes. stops shared_ptr from initializing to nullptr. 
    Mesh(std::size_t dims) : m_size(dims){} 

    // forward args to std::vector
    template<typename ArgType>
    Mesh(const Eigen::MatrixBase<ArgType>& xpr, std::size_t dims=1)
      : m_size(dims)
    {
      if(dims > numDimsMax) throw std::runtime_error("Can't construct mesh with given dims"); 
      for(std::size_t idx=0; idx<m_size; ++idx){
        m_mesh_arr[idx] = xpr; 
      }
    }

    template<typename... ArgType>
    Mesh(const Eigen::MatrixBase<ArgType>&... xpr)
      : m_size(0)
    {
      static_assert(sizeof...(ArgType), "Can't construct mesh with given dims");
      auto lam = [&](const auto& x){m_mesh_arr[m_size] = x; ++m_size; }; 
      std::apply(
        lam, 
        xpr... 
      ); 
    }

    // Copy 
    Mesh(const Mesh& other)=default; 
    
    // destructors
    ~Mesh()=default;

    // Member Functions ===============================================================

    // get list of stored 1d mesh pointers. uses NullaryExpr + Converter to propragate constness 
    auto& getAxesList(){ return m_mesh_arr; } 
    
    // number of dimensions 
    std::size_t numDims() const {return m_size; } 

    // size of a specific axis 
    std::size_t sizeOfDim(std::size_t i) const {return m_mesh_arr[i].size();} 

    // get a specific axis
    auto& getAxis(std::size_t i){ return (m_mesh_arr[i]); }
    const auto& getAxis(std::size_t i) const { return (m_mesh_arr[i]); }

    // get a spec
    auto& getAxisSafe(std::size_t i)
    {
      if(i >= m_size) throw std::runtime_error("error ith_dim out of range. ");  
      return ( m_mesh_arr.at(i) ); 
    }
    const auto& getAxisSafe(std::size_t i) const 
    {
      if(i >= m_size) throw std::runtime_error("error ith_dim out of range. ");  
      return ( m_mesh_arr.at(i) ); 
    }

    // full size of XD mesh. i.e. axis1.size() * ... * axisn.size()
    std::size_t sizesProduct() const 
    {
      std::size_t p = 1; 
      for(auto idx=0; idx<m_size; ++idx){
        p *= m_mesh_arr[idx].size(); 
      }
      return p; 
    } 

    // product of axes up to dim exclusively [first, dim)
    std::size_t sizesMiddleProduct(std::size_t start, std::size_t end) const 
    {
      if(start > end) throw std::invalid_argument("start index must be <= end index for middle product"); 
      if(end > m_size) throw std::invalid_argument("end index must be <= # of numDims in MeshXD"); 
      std::size_t prod = 1; 
      for(auto idx=start; idx<end; ++idx){
        prod *= m_mesh_arr[idx].size(); 
      }
      return prod; 
    } 

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

  private:
    // Nested FO Structs
    struct Converter
    {
      const std::vector<std::shared_ptr<Eigen::VectorXd>>& m_wrapped; 
      const Eigen::VectorXd* operator()(Eigen::Index i) const { return m_wrapped[i].get(); }
    }; 
};

template<typename... Args> 
auto make_Mesh(Args... args)
{
  return std::make_shared<Mesh>(args...); 
}

} // end namespace fdm 

#endif // Mesh.hpp 