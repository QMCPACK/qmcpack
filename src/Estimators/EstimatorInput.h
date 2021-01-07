//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESIMATORINPUT_H
#define QMCPLUSPLUS_ESIMATORINPUT_H

#include <variant>

#include "SpinDensityInput.h"

namespace qmcplusplus
{
/** Will be an variant over specific estimators input structures
 *
 *  Why not inheritance? Because this is a closed set which each input having the
 *  interface its Estimator needs for its input. With separate input parsing and 
 *  simulation object costruction phases the estimator manager will need the list of
 *  input structures to construct its estimators, this can be done in a type safe manner
 *  using a vector of EstimatorInputs.  
 *
 *  The match between the input estimator type string and estimator input type will only happen in the input parsing.
 *  
 *  using EstimatorInput = std::variant<SpinDensityInput,DensityMatrices1BInput,...>;
 */
using EstimatorInput = SpinDensityInput;

} // namespace qmcplusplus

#endif
