//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
//////////////////////////////////////////////////////////////////////////////////////

#include "WaveFunctionTesterBatchedInput.h"

#include "Message/UniformCommunicateError.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
void WaveFunctionTesterBatchedInput::readXML(xmlNodePtr xml_input)
{
  std::string basic{"yes"};
  std::string ratio{"no"};

  ParameterSet parameter_set;
  parameter_set.add(basic, "basic", {"yes", "no"});
  parameter_set.add(ratio, "ratio", {"yes", "no"});
  parameter_set.add(delta_, "delta");
  parameter_set.add(tolerance_, "tolerance");
  parameter_set.put(xml_input);

  for (xmlNodePtr child = xml_input->children; child != nullptr; child = child->next)
  {
    if (std::string(castXMLCharToChar(child->name)) != "parameter")
      continue;

    const std::string parameter_name(getXMLAttributeValue(child, "name"));
    if (parameter_name == "clone" || parameter_name == "hamiltonianpbyp" || parameter_name == "source" ||
        parameter_name == "orbitalutility" || parameter_name == "printEloc" || parameter_name == "virtual_move" ||
        parameter_name == "sd")
      throw UniformCommunicateError("Batched wavefunction tester does not support parameter \"" + parameter_name +
                                    "\".");
  }

  run_basic_ = basic == "yes";
  run_ratio_ = ratio == "yes";
  if (!run_basic_ && !run_ratio_)
    throw UniformCommunicateError("Batched wavefunction tester requires basic=\"yes\" or ratio=\"yes\".");
}
} // namespace qmcplusplus