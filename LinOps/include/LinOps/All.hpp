// All.hpp
// 
// header file to include full linop framework
// 
// JAF 12/5/2025 

#ifndef LINOP_ALL_H
#define LINOP_ALL_H

#include "Mesh1D.hpp"
#include "Vector1D.hpp"
#include "MeshXD.hpp"
#include "VectorXD.hpp"

#include "LinOpTraits.hpp"
#include "LinOpMixIn.hpp"
#include "LinOpBase.hpp"
#include "LinOpExpr.hpp"

#include "Operators/IOp.hpp"
#include "Operators/RandOp.hpp"
#include "Operators/DirectionalRandOp.hpp"

#include "DiffOps/NthDerivOp.hpp"
#include "DiffOps/DirectionalNthDerivOp.hpp" 

#include "CoeffMixIn.hpp"
#include "Coeffs/AutonomousCoeff.hpp"
#include "Coeffs/TimeDepCoeff.hpp"

#endif // All.hpp
