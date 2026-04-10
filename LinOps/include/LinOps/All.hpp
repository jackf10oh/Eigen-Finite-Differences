// All.hpp
// 
// header file to include full linop framework
// 
// JAF 12/5/2025 

#ifndef LINOP_ALL_H
#define LINOP_ALL_H

#include "Mesh1D.hpp"
#include "Vector1D.hpp"
#include "VectorXD.hpp"
#include "MeshXD.hpp"

#include "LinOpTraits.hpp"
#include "LinOpBase.hpp"
#include "LinOpExpr.hpp"

#include "Operators/IOp.hpp"
#include "Operators/RandOp.hpp"
#include "Operators/DirectionalRandOp.hpp"
#include "Operators/DiffOps/NthDerivOp.hpp"
#include "Operators/DiffOps/DirectionalNthDerivOp.hpp" 

#include "CoeffOpMixIn.hpp"
#include "Operators/CoeffOps/AutonomousCoeff.hpp"
#include "Operators/CoeffOps/TimeDepCoeff.hpp"

#endif // All.hpp
