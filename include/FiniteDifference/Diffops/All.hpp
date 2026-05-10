// All.hpp 
//
// include everything from Diffops 
// 
// JAF 5/3/2026 

#ifndef FDM_DIFFOPS_ALL_H
#define FDM_DIFFOPS_ALL_H

#ifndef EIGEN_SPARSEMATRIXBASE_PLUGIN
  #define EIGEN_SPARSEMATRIXBASE_PLUGIN <FiniteDifference/EigenFdmPlugin.hpp> 
#endif 

#include "Traits.hpp"
#include "Functors.hpp"
#include "PartialDerivBase.hpp"
#include "EigenEvaluator.hpp"
#include "EvaluatorBase.hpp"
#include "NwiseBinaryOp.hpp" 
#include "NwiseUnaryOp.hpp" 
#include "NthPartialDeriv.hpp" 

#endif // All.hpp (Diffops)