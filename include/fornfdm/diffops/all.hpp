// all.hpp 
//
// include everything from diffops 
// 
// JAF 5/3/2026 

#ifndef FORNFDM_DIFFOPS_ALL_H
#define FORNFDM_DIFFOPS_ALL_H

#ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
  #define EIGEN_SPARSEMATRIXBASE_PLUGIN <FornFdm/EigenFdmPlugin.hpp> 
#endif 

#include "traits.hpp"
#include "Functors.hpp"
#include "PartialDerivBase.hpp"
#include "KroneckerEvaluator.hpp"
#include "EvaluatorBase.hpp"
#include "NwiseBinaryOp.hpp" 
#include "NwiseUnaryOp.hpp" 
#include "NthPartialDeriv.hpp" 

#endif // all.hpp (diffops)