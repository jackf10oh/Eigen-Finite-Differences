// all.hpp 
//
// Includs all files from every subdirectory,
// except for utilitis. 
//
// JAF 4/11/2026 

#ifndef FORNFDM_ALL_H
#define FORNFDM_ALL_H

#ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
  #define EIGEN_SPARSEMATRIXBASE_PLUGIN <FornFdm/EigenFdmPlugin.hpp> 
#endif 

#include "Types.hpp"
#include "Traits.hpp"
#include "Mesh.hpp" 

#include "diffops/all.hpp"
#include "coeffs/all.hpp"
#include "outside_steps/all.hpp"
#include "texprs/all.hpp"
#include "Solvers/all.hpp"

#endif

