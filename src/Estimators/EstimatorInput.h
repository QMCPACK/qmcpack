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
/** Input variant for Estimators
 *
 *  Why not inheritance? Because this is a closed set which each input having the
 *  interface its Estimator needs for its input.
 */
  using EstimatorInput = std::variant<SpinDensityInput>;

}

#endif
