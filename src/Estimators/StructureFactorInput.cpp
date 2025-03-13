//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SkEstimator, SkAllEstimator
//////////////////////////////////////////////////////////////////////////////////////

#include "StructureFactorInput.h"

namespace qmcplusplus
{

StructureFactorInput::StructureFactorInput(xmlNodePtr cur)
{
  // This results in checkParticularValidity being called on StructureFactorInput
  input_section_.readXML(cur);

  auto setIfInInput = [&](auto& var, const std::string& tag) -> bool { return input_section_.setIfInInput(var, tag); };
  setIfInInput(name_, "name");
  setIfInInput(type_, "type");
  setIfInInput(write_hdf5_, "writehdf5");
  setIfInInput(write_rho_, "writerho");
  setIfInInput(write_ion_ion_, "writeionion");
  setIfInInput(source_, "source");
  setIfInInput(target_, "target");
}

void StructureFactorInput::StructureFactorInputSection::checkParticularValidity()
{
  const std::string error_tag{"StructureFactor input: "};
  if(has("writerho"))
    if(get<bool>("writerho"))
      if(!(has("source") && has("target")))
	throw UniformCommunicateError(error_tag + " writerho requires explicit source and target definition!");

  if(has("writeionion"))
    if(get<bool>("writeionion"))
      if(!(has("source") && has("target")))
	throw UniformCommunicateError(error_tag + " writerionion requires explicit source and target definition!");
}

}
