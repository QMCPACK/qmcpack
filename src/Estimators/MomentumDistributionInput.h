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
class MomentumDistributionInput : public InputSection
{
public:
  MomentumDistributionInput()
  {
    section_name   = "MomentumDistribution";
    attributes     = {"type","name","samples","kmax"};
    strings        = {"type","name"};
    integers       = {"samples"};
    reals          = {"kmax"};
    default_values = { {"samples", 40} };
  };
};

} // namespace qmcplusplus
#endif /* MOMENTUMDISTRIBUTIONINPUT_H */
