// all.hpp 
//
// Includs all files from every subdirectory,
// except for utilitis. 
//
// JAF 4/11/2026 

#ifndef FORNFDM_ALL_H
#define FORNFDM_ALL_H

#ifndef FORNFDM_PLUGIN_SET
  #define FORNFDM_PLUGIN_SET
  #ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
    #define EIGEN_SPARSEMATRIXBASE_PLUGIN <fornfdm/plugin.hpp> 
  #else
    #error "fornfdm requires usage of EIGEN_SPARSEMATRIXBASE_PLUGIN macro. move other macros to EIGEN_SPARSEMATRIXBASE_PLUGIN_OTHER"
  #endif
#endif

#include "types.hpp"
#include "traits.hpp"
#include "Mesh.hpp" 

#include "diffops/all.hpp"
#include "coeffs/all.hpp"
#include "outside_steps/all.hpp"
#include "texprs/all.hpp"
#include "solvers/all.hpp"

#endif

