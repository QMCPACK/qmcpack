//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from:  QMCHamiltonians/SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "SpaceGridInput.h"

#include <sstream>
#include <vector>

#include "EstimatorInput.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

SpaceGridInput::SpaceGridAxisInput::SpaceGridAxisInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(label_, "label");
  setIfInInput(grid_, "grid");
  setIfInInput(p1_, "p1");
  setIfInInput(p2_, "p2");
  setIfInInput(scale_, "scale");
}

void SpaceGridInput::SpaceGridAxisInput::SpaceGridAxisInputSection::setFromStreamCustom(const std::string& ename,
                                                                                        const std::string& name,
                                                                                        std::istringstream& svalue)
{
  // grid can only be an attribute for space grid axes.
  if (name == "grid")
  {
    try
    {
      values_[name] = parseGridInput<Real>(svalue);
    }
    catch (const UniformCommunicateError& uce)
    {
      std::ostringstream msg;
      msg << "SpaceGridAxisInputSection failed in custom stream handler with: " << uce.what()
          << " a report of the current parse progress should be found above\n";
      report();
      throw UniformCommunicateError(msg.str());
    }
  }
  else
    throw std::runtime_error("bad name passed or custom setFromStream not implemented in derived class.");
}

SpaceGridInput::SpaceGridOriginInput::SpaceGridOriginInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(p1_, "p1");
  setIfInInput(p2_, "p2");
  setIfInInput(fraction_, "fraction");
}

SpaceGridInput::SpaceGridInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(coord_form_, "coord");
  // rip open the axes inputs we're guarateed they have the proper dimensionality already

  auto axes = input_section_.get<std::vector<std::any>>("axis");
  checkAxes(axes);
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    auto* axis_input = std::any_cast<SpaceGridAxisInput>(&axes[d]);
    axis_labels_[d]  = axis_input->get_label();
    axis_p1s_[d]     = axis_input->get_p1();
    axis_p2s_[d]     = axis_input->get_p2();
    axis_grids_[d]   = axis_input->get_grid();
    axis_scales_[d]  = axis_input->get_scale();
  }
  checkGrids();
  if (input_section_.has("origin"))
  {
    auto space_grid_origin = input_section_.get<SpaceGridOriginInput>("origin");
    origin_p1_             = space_grid_origin.get_p1();
    origin_p2_             = space_grid_origin.get_p2();
    origin_fraction_       = space_grid_origin.get_fraction();
  }
}

void SpaceGridInput::checkAxes(std::vector<std::any>& axes)
{
  auto& ax_labels = axes_label_sets.at(coord_form_);
  for (auto& axis : axes)
  {
    auto* axis_input       = std::any_cast<SpaceGridAxisInput>(&axis);
    std::string axis_label = axis_input->get_input().template get<std::string>("label");
    auto result            = std::find(std::begin(ax_labels), std::end(ax_labels), axis_label);
    if (result == std::end(ax_labels))
      throw UniformCommunicateError(axis_label + " is not a valid label for coord form " +
                                    input_section_.get<std::string>("coord"));
  }
}

void SpaceGridInput::checkGrids()
{
  //check that all axis grid values fall in the allowed intervals for the coord label
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    if (axis_labels_[d] == "r" || axis_labels_[d] == "phi" || axis_labels_[d] == "theta")
    {
      if (axis_grids_[d].umin < 0.0 || axis_grids_[d].umax > 1.0)
      {
        std::ostringstream error;
        error << " grid values for " << axis_labels_[d] << " must fall in [0,1]" << std::endl;
        error << " interval provided: [" << axis_grids_[d].umin << "," << axis_grids_[d].umax << "]" << std::endl;
        throw UniformCommunicateError(error.str());
      }
    }
    // all other legal labels {"x","y","z"} can be over -1.0 to 1.0 
    else if (axis_grids_[d].umin < -1.0 || axis_grids_[d].umax > 1.0)
    {
      std::ostringstream error;
      error << "  grid values for " << axis_labels_[d] << " must fall in [-1,1]" << std::endl;
      error << "  interval provided: [" << axis_grids_[d].umin << "," << axis_grids_[d].umax << "]" << std::endl;
      throw UniformCommunicateError(error.str());
    }
  }
}


std::any SpaceGridInput::SpaceGridInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}

std::any makeSpaceGridInput(xmlNodePtr cur, std::string& value_label)
{
  SpaceGridInput sgi{cur};
  value_label = "spacegrid";
  return sgi;
}

void SpaceGridInput::SpaceGridInputSection::checkParticularValidity()
{
  auto axes = std::any_cast<std::vector<std::any>>(values_["axis"]);
  static_assert(std::is_same<decltype(axes), std::vector<std::any>>::value);
  auto axis_count = axes.size();
  if (axis_count != OHMMS_DIM)
  {
    std::ostringstream error;
    error << "SpaceGrid input must contain " << OHMMS_DIM << " axes, " << axis_count << " found!";
    throw UniformCommunicateError(error.str());
  }
}

} // namespace qmcplusplus
