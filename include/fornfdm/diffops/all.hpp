// all.hpp 
//
// include everything from diffops 
// 
// JAF 5/3/2026 

#ifndef FORNFDM_DIFFOPS_ALL_H
#define FORNFDM_DIFFOPS_ALL_H

// #include "../plugin.hpp" // not forced in includes
#include "traits.hpp"
#include "functors.hpp"
#include "NodeSelector.hpp"
#include "CenteredNodeSelector.hpp"
#include "ForwardNodeSelector.hpp"
#include "BackwardNodeSelector.hpp"
#include "PeriodicNodeSelector.hpp"
#include "EvaluatorBase.hpp"
#include "TimeEvaluation.hpp"
#include "PartialDerivBase.hpp"
#include "EigenEvaluatorImpl.hpp"
#include "NthPartialDeriv.hpp"
#include "NwiseUnaryOp.hpp"
#include "NwiseBinaryOp.hpp"
#include "MatrixFree.hpp"

#endif // all.hpp (diffops)