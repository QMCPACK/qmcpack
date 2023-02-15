//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////
#include "MagnetizationDensityInput.h"

#include <cmath>

#include "OhmmsData/AttributeSet.h"
#include "Message/UniformCommunicateError.h"
#include "EstimatorInput.h"
namespace qmcplusplus
{

MagnetizationDensityInput::MagnetizationDensityInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(nsamples_, "samples");
  setIfInInput(integrator_,"integrator"); 
}


MagnetizationDensityInput::DerivedParameters MagnetizationDensityInput::calculateDerivedParameters(const Lattice& lattice) const
{
  PosType corner = 0.0;
  if (have_center_)
    corner = center_ - lattice.Center;
  else if (have_corner_)
    corner = corner_;

  TinyVector<int, DIM> grid;
  if (have_dr_)
    for (int d = 0; d < DIM; ++d)
      grid[d] = (int)std::ceil(std::sqrt(dot(lattice.Rv[d], lattice.Rv[d])) / dr_[d]);
  else if (have_grid_)
    grid = grid_;

  size_t npoints = 1;
  for (int d = 0; d < DIM; ++d)
    npoints *= grid[d];

  TinyVector<int, DIM> gdims;
  gdims[0] = npoints / grid[0];
  for (int d = 1; d < DIM; ++d)
    gdims[d] = gdims[d - 1] / grid[d];

  return {corner, grid, gdims, npoints};
}

std::any MagnetizationDensityInput::MagnetizationDensityInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}

}//namespace qmcplusplus
