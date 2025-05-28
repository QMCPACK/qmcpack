//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EnergyDensityInput.h"
#include "EstimatorInput.h"

namespace qmcplusplus
{

EnergyDensityInput::EnergyDensityInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(name_, "name");
  setIfInInput(type_, "type");
  setIfInInput(dynamic_, "dynamic");
  setIfInInput(static_, "static");
  setIfInInput(ion_points_, "ion_points");
  setIfInInput(ref_points_input_, "reference_points");
}

std::vector<SpaceGridInput> EnergyDensityInput::get_space_grid_inputs() const
{
  auto any_vec = input_section_.get<std::vector<std::any>>("spacegrid");
  std::vector<SpaceGridInput> sgis;
  for (auto& grid : any_vec)
    sgis.emplace_back(std::any_cast<SpaceGridInput>(grid));
  return sgis;
}

} // namespace qmcplusplus
