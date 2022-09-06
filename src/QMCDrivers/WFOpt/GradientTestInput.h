//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_GRADIENTTESTINPUT_H
#define QMCPLUSPLUS_GRADIENTTESTINPUT_H

#include "io/OhmmsData/libxmldefs.h"

namespace qmcplusplus
{

class GradientTestInput
{
public:
  GradientTestInput() {}
  void readXML(xmlNodePtr xml_input);

protected:
  bool do_param_output_ = false;
  double finite_diff_delta_ = 1e-5;

public:
  bool do_param_output() const { return do_param_output_; }
  double get_finite_diff_delta() const { return finite_diff_delta_; }
};

} // namespace qmcplusplus

#endif
