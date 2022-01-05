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

#include "GradientTestInput.h"
#include "io/OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
void GradientTestInput::readXML(xmlNodePtr xml_input)
{
  ParameterSet param;
  std::string do_output_file = "no";
  param.add(do_output_file, "output_param_file");
  param.put(xml_input);

  do_param_output_ = (do_output_file != "no");
}


} // namespace qmcplusplus
