// NodeSelector.hpp 
//
// Class responsible for selecting finite-difference scheme
// i.e. forward stencil, centered stencil, backward stencil 
//
// JAF 5/10/2026 

#ifndef FDM_DIFFOPS_NODESELECTOR_H
#define FDM_DIFFOPS_NODESELECTOR_H

namespace fdm{
namespace linops{
namespace internal{

// holds (x1, x2, ..., xn) in the same dimension 
template<class tag, std::size_t numNodesMin>
struct NodeSelector{}; 

// TODO struct Forward {}; // struct Backward {};

} // end namespace internal 
} // end namespace linops 
} // end namespace fdm

#include "CenteredNodeSelector.hpp"
#include "ForwardNodeSelector.hpp"
#include "BackwardNodeSelector.hpp"

#endif // NodeSelector.hpp 