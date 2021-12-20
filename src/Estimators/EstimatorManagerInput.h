//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESIMATORMANAGERINPUT_H
#define QMCPLUSPLUS_ESIMATORMANAGERINPUT_H

#include "Configuration.h"
#include "InputSection.h"

namespace qmcplusplus
{
/** Input representation for Driver base class runtime parameters
 */
class EstimatorManagerInput
{
public:
  class EstimatorMangerInputSection : public InputSection
  {
  public:
    /** parse time definition of input parameters */
    EstimatorManagerInputSection()
    {
      // clang-format off
      section_name  = "OneBodyDensityMatrix";
      attributes    = {"name", "type"};
      parameters    = {"basis", "energy_matrix", "integrator", "evaluator", "scale",
                       "corner", "center", "points", "samples", "warmup", "timestep",
                       "use_drift", "check_overlap", "check_derivatives", "acceptance_ratio", "rstats",
                       "normalized", "volumed_normed"};
      bools         = {"energy_matrix", "use_drift", "normalized", "volume_normed",
                       "check_overlap", "check_derivatives", "rstats", "acceptance_ratio"};
      enums         = {"integrator", "evaluator"};
      strings       = {"name", "type"};
      multi_strings = {"basis"};
      integers      = {"points", "samples"};
      reals         = {"scale", "timestep"};
      positions     = {"center", "corner"};
      required      = {"name", "basis"};
      // I'd much rather see the default defined in simple native c++ as below
      // clang-format on
    }
  };
};
}
