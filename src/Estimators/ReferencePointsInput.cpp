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
// Some code refactored from:  QMCHamiltonians/ReferencePoints.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "ReferencePointsInput.h"

#include <sstream>

#include "EstimatorInput.h"
#include "ModernStringUtils.hpp"

namespace qmcplusplus
{
template bool InputSection::setIfInInput<ReferencePointsInput::Coord>(ReferencePointsInput::Coord& var,
                                                                      const std::string& tag);

ReferencePointsInput::ReferencePointsInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(coord_form_, "coord");
  readRefPointsXML(cur);
}

// ReferencePointsInput::ReferencePointsInput(const Points& points, const CoordForm coord_form)
//     : points_(points), coord_form_(coord_form)
// {}

void ReferencePointsInput::readRefPointsXML(xmlNodePtr cur)
{
  using modernstrutil::split;
  using modernstrutil::strip;

  // read refpoints values they are some sequence of value nodes
  std::string node_str{XMLNodeString{cur}};
  std::vector<std::string_view> lines = split(strip(node_str), "\n");
  for (int i = 0; i < lines.size(); i++)
  {
    auto stripped_line                   = strip(lines[i]);
    std::vector<std::string_view> tokens = split(stripped_line, " ");
    if (tokens.size() != OHMMS_DIM + 1)
    {
      std::ostringstream error;
      error << error_tag << "reference point has 4 entries, given " << tokens.size() << ": " << lines[i];
      throw UniformCommunicateError(error.str());
    }
    else
    {
      Point rp;
      for (int d = 0; d < OHMMS_DIM; d++)
      {
        try
        {
          rp[d] = std::stod(std::string(tokens[d + 1].begin(), tokens[d + 1].size()));
        }
        catch (const std::invalid_argument& ia)
        {
          throw UniformCommunicateError(ia.what());
        }
      }

      // This must be done in constructor of ReferencePoints
      // rp                = dot(crd, rp);
      points_[std::string(tokens[0].begin(), tokens[0].size())] = rp;
    }
  }
}

std::any ReferencePointsInput::ReferencePointsInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}

std::any makeReferencePointsInput(xmlNodePtr cur, std::string& value_label)
{
  ReferencePointsInput rpi{cur};
  value_label = "referencepoints";
  return rpi;
}

} // namespace qmcplusplus
