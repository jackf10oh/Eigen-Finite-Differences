// All.hpp 
//
// Includs all files from every subdirectory,
// except for utilitis. 
//
// JAF 4/11/2026 

#ifndef FORNFDM_ALL_H
#define FORNFDM_ALL_H

#ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
  #define EIGEN_SPARSEMATRIXBASE_PLUGIN <FiniteDifference/EigenFdmPlugin.hpp> 
#endif 

#include "Types.hpp"
#include "Traits.hpp"
#include "Mesh.hpp" 

#include "Diffops/All.hpp"
#include "Coeffs/All.hpp"
#include "OutsideSteps/All.hpp"
#include "TExprs/All.hpp"
#include "Solvers/All.hpp"

#endif

