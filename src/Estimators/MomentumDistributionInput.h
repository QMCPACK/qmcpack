//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: MomentumEstimator.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_MOMENTUMDISTRIBUTIONINPUT_H
#define QMCPLUSPLUS_MOMENTUMDISTRIBUTIONINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{


/** Native representation for Momentum Distribution Estimators inputs
 */
// clang-format: off
class MomentumDistributionInput : public InputSection
{
public:
  MomentumDistributionInput()
  {
    section_name = "MomentumDistribution";
    attributes   = {"type", "name", "samples", "kmax", "kmax0", "kmax1", "kmax2"};
    strings      = {"type", "name"};
    integers     = {"samples"};
    reals        = {"kmax", "kmax0", "kmax1", "kmax2"};
    default_values =
        {{"name"   , std::string("nofk")}, 
         {"samples", int(40)}, 
         {"kmax"   , Real(0.0)}, 
         {"kmax0"  , Real(0.0)}, 
         {"kmax1"  , Real(0.0)}, 
         {"kmax2"  , Real(0.0)}};
  };
};
// clang-format: on

} // namespace qmcplusplus
#endif /* MOMENTUMDISTRIBUTIONINPUT_H */
