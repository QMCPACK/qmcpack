//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EnergyDensityInput.h"

#include "string_utils.h"
#include "EstimatorInput.h"

namespace qmcplusplus
{

EnergyDensityInput::EnergyDensityInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(name_, "name");
  // setIfInInput(type_, "type");
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


// bool EnergyDensityInput::EnergyDensityInputSection::setFromStreamCustom(const std::string& name, istringstream sstream)
// {
//   ReferencePoints& rp = ref_points_;
//   bool succeeded = true;
//   put(P, Pref);
//   OhmmsAttributeSet ra;
//   std::string coord = "";
//   ra.add(coord, "coord");
//   ra.put(cur);
//   for (int i = 0; i < DIM; i++)
//     for (int d = 0; d < DIM; d++)
//       axes(d, i) = P.getLattice().a(i)[d];
//   Tensor_t crd;
//   if (coord == "cell")
//   {
//     coordinate = cellC;
//     crd        = axes;
//   }
//   else if (coord == "cartesian")
//   {
//     coordinate = cartesianC;
//     for (int i = 0; i < DIM; i++)
//       for (int d = 0; d < DIM; d++)
//         if (d == i)
//           crd(i, i) = 1.0;
//         else
//           crd(d, i) = 0.0;
//   }
//   else
//   {
//     app_log() << std::endl;
//     app_log() << "    Valid coordinates must be provided for element reference_points." << std::endl;
//     app_log() << "      You provided: " << coord << std::endl;
//     app_log() << "      Options are cell or cartesian." << std::endl;
//     app_log() << std::endl;
//     succeeded = false;
//   }
//   //read in the point contents
//   app_log() << "    reading reference_points contents" << std::endl;
//   std::vector<std::string> lines = split(strip(XMLNodeString{cur}), "\n");
//   for (int i = 0; i < lines.size(); i++)
//   {
//     std::vector<std::string> tokens = split(strip(lines[i]), " ");
//     if (tokens.size() != DIM + 1)
//     {
//       app_log() << "  reference point has 4 entries, given " << tokens.size() << ": " << lines[i] << std::endl;
//       succeeded = false;
//     }
//     else
//     {
//       Point rp;
//       for (int d = 0; d < DIM; d++)
//       {
//         rp[d] = string2real(tokens[d + 1]);
//       }
//       rp                = dot(crd, rp);
//       points[tokens[0]] = rp;
//     }
//   }
//   return succeeded;
// }

} // namespace qmcplusplus
