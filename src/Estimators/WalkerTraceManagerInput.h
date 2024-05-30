//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_WALKERTRACEMANAGERINPUT_H
#define QMCPLUSPLUS_WALKERTRACEMANAGERINPUT_H

#include "InputSection.h"

namespace qmcplusplus
{

/** Native representation for WalkerTraceManager input
 */
struct WalkerTraceInput : public InputSection
{
  bool present;

  WalkerTraceInput() : WalkerTraceInput(NULL) {}

  WalkerTraceInput(xmlNodePtr cur) 
  {
    section_name   = "walkertraces";
    attributes     = {"step_period", "particle", "min", "max", "median", "qtiles", "verbose"};
    integers       = {"step_period"};
    bools          = {"particle", "min", "max", "median", "qtiles", "verbose"};
    default_values = {{"step_period", int(1)},
                      {"particle", bool(false)},
                      {"min"     , bool(true)},
                      {"max"     , bool(true)},
                      {"median"  , bool(true)},
                      {"qtiles"  , bool(true)},
                      {"verbose" , bool(false)}};
    present = cur!=NULL;
    if(present)
      readXML(cur);
  };
};

} // namespace qmcplusplus
#endif /* WALKERTRACEMANAGERINPUT_H */
