//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_WALKERLOGINPUT_H
#define QMCPLUSPLUS_WALKERLOGINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

/** Native representation for walker logs input
 */
struct WalkerLogInput : public InputSection
{
  bool present;

  WalkerLogInput() : WalkerLogInput(NULL) {}

  WalkerLogInput(xmlNodePtr cur)
  {
    section_name   = "walkerlogs";
    attributes     = {"step_period", "particle", "min", "max", "median", "quantiles", "verbose"};
    integers       = {"step_period"};
    bools          = {"particle", "min", "max", "median", "quantiles", "verbose"};
    default_values = {{"step_period", int(1)}, {"particle", bool(false)}, {"min", bool(true)},     {"max", bool(true)},
                      {"median", bool(true)},  {"quantiles", bool(true)}, {"verbose", bool(false)}};
    present        = cur != NULL;
    if (present)
      readXML(cur);
  };
};

} // namespace qmcplusplus
#endif /* WALKERLOGINPUT_H */
