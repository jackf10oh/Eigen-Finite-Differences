// all.hpp 
//
// include everything from diffops 
// 
// JAF 5/3/2026 

#ifndef FORNFDM_DIFFOPS_ALL_H
#define FORNFDM_DIFFOPS_ALL_H

#ifndef FORNFDM_PLUGIN_SET
  #define FORNFDM_PLUGIN_SET
  #ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
    #define EIGEN_SPARSEMATRIXBASE_PLUGIN <fornfdm/plugin.hpp> 
  #else
    #error "fornfdm requires usage of EIGEN_SPARSEMATRIXBASE_PLUGIN macro. move other macros to EIGEN_SPARSEMATRIXBASE_PLUGIN_OTHER"
  #endif
#endif

#include "traits.hpp"
#include "functors.hpp"
#include "PartialDerivBase.hpp"
#include "KroneckerEvaluator.hpp"
#include "EvaluatorBase.hpp"
#include "NwiseBinaryOp.hpp" 
#include "NwiseUnaryOp.hpp" 
#include "NthPartialDeriv.hpp" 

#endif // all.hpp (diffops)