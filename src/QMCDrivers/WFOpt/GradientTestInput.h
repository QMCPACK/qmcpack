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

public:
  bool do_param_output() { return do_param_output_; }
};

} // namespace qmcplusplus
