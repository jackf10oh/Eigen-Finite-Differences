// NthDerivOp.hpp
//
//
//
// JAF 12/5/2025

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<cstdint>
#include "../../Utilities/FornbergCalc.hpp"
#include "../LinOpMixIn.hpp"
#include "../LinOpBase.hpp" 

namespace fdm{
  namespace linops{

template<std::size_t orderN>
class NthDerivOp : public LinOpMixIn<NthDerivOp<orderN>>, public LinOpBase1D<NthDerivOp<orderN>>
{
  private:
    // member data ----------------------------------------------------------- 
    fdm::Matrix m_stencil; 
  public:
    static constexpr std::size_t order = orderN; 

    // Constructors + Destructor ===================================================
    // default 
    NthDerivOp()=default; 

    // from mesh
    NthDerivOp(const fdm::SharedConstMesh1D& m)
    {this->setMesh(m);};

    // copy 
    NthDerivOp(const NthDerivOp& other)=default; 

    // destructor
    ~NthDerivOp()=default; 
    
    // Member Funcs =====================================================

    // Matrix Getters 
    const fdm::Matrix& asMatrix() const { return m_stencil; };  
    
  protected: 
    // Unreachable ------------------------------------------------------------
    // set the mesh domain the derivative operator works on 
    void setMesh1D_impl(const fdm::SharedConstMesh1D& m)
    {
      if constexpr(order == 1){ 
        m_stencil.setIdentity(); 
        return; 
      } 

      const std::size_t mesh_size = m->size();

      // resize matrix to fit
      m_stencil.resize(mesh_size,mesh_size);
      
      constexpr int one_sided_skirt = order;  
      constexpr int centered_skirt = (order+1)/2;  

      // allocate full list of coeff triples 
      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.resize(2*centered_skirt*(1+orderN) + (mesh_size-2*centered_skirt)*(1+2*centered_skirt)); // not sure if this is correct size of lise...

      // begin OpenMP parallel section. **** COMMENTED OUT: std::vector writing is not thread safe ******
      // #pragma omp parallel 
      {
        // instantiate stateful fornberg calculator. one per thread 
        fdm::utils::FornCalc weight_calc(1+2*((orderN+1)/2),orderN);
        // first rows with forward stencil
        // #pragma omp for nowait 
        for(std::size_t i=0; i<centered_skirt; i++)
        {
          auto left = m->cbegin()+i;
          auto right = left+one_sided_skirt+1; 
          auto weights = weight_calc.GetWeights(*left, left,right, orderN); 
          std::size_t offset=i; 
          for(auto& w : weights){
            tripletList[i*(1+one_sided_skirt) 
                        + offset] = T(i,offset,w);
            offset++; 
          }
        }
        // middle rows with centered centered stencil  
        // #pragma omp for nowait 
        for(std::size_t i=centered_skirt;i<mesh_size-centered_skirt; i++)
        {
          auto left = m->cbegin()-centered_skirt+i;  
          auto right = left+(centered_skirt+1+centered_skirt);
          auto weights = weight_calc.GetWeights(m->at(i), left,right, orderN);
          int offset = -centered_skirt;
          for(auto& w : weights){
            tripletList[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = T(i,i+offset,w);
            offset++; 
          };
        }
        // last rows 
        // #pragma omp for nowait 
        for(std::size_t i = mesh_size-centered_skirt; i<mesh_size; i++)
        {
          auto right = m->cbegin() + i; 
          auto left = right-(one_sided_skirt+1); 
          auto weights = weight_calc.GetWeights(*right, left,right, orderN); 
          int offset= -one_sided_skirt;
          for(auto it=weights.cbegin(); it!=weights.cend(); it++){
            tripletList[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = T(i,i+offset,*it);  
            offset++;
          }
        }
      } 
      // End OpenMP parallel section. implicit barrier 
      
      // Eigen 5.0 :(. might be faster if sorting time > time saved from parallelism 
      // m_stencil.insertFromSortedTriplets()(tripletList.begin(), tripletList.end());   
      m_stencil.setFromTriplets(tripletList.begin(), tripletList.end()); 
    }

}; 

  } // end namespace linops 
} // end namespace fdm  
#endif // NthDerivOp.hpp
