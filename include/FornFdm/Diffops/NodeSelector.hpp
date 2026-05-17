// NodeSelector.hpp 
//
// Class responsible for selecting finite-difference scheme
// i.e. forward stencil, centered stencil, backward stencil 
//
// JAF 5/10/2026 

#ifndef FORNFDM_DIFFOPS_NODESELECTOR_H
#define FORNFDM_DIFFOPS_NODESELECTOR_H

namespace fornfdm{
namespace linops{
namespace internal{

// holds (x1, x2, ..., xn) in the same dimension 
template<class tag, std::size_t numNodesMin>
struct NodeSelector{}; 

} // end namespace internal 
} // end namespace linops 
} // end namespace fornfdm

#include "CenteredNodeSelector.hpp"
#include "ForwardNodeSelector.hpp"
#include "BackwardNodeSelector.hpp"

#endif // NodeSelector.hpp 