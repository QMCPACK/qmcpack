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
  setIfInInput(integrator_, "integrator");
  have_center_ = setIfInInput(center_, "center");
  have_corner_ = setIfInInput(corner_, "corner");
  have_grid_   = setIfInInput(grid_real_, "grid");
  have_dr_     = setIfInInput(dr_, "dr");
}


MagnetizationDensityInput::DerivedParameters MagnetizationDensityInput::calculateDerivedParameters(
    const Lattice& lattice) const
{
  PosType corner = 0.0;
  //Corner and center can be taken care of by defaults.  corner=0 if not specified.
  if (have_center_)
    corner = center_ - lattice.Center;
  else if (have_corner_)
    corner = corner_;

  TinyVector<int, DIM> grid;

  //dr or grid must be specified to perform the grid math.  Input already checked before we get here.
  if (have_dr_)
    for (int d = 0; d < DIM; ++d)
      grid[d] = (int)std::ceil(std::sqrt(dot(lattice.Rv[d], lattice.Rv[d])) / dr_[d]);
  else if (have_grid_)
    for (int d = 0; d < DIM; ++d)
      grid[d] = (int)std::ceil(grid_real_[d]);

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

void MagnetizationDensityInput::MagnetizationDensityInputSection::checkParticularValidity()
{
  using namespace estimatorinput;
  const std::string error_tag{"MagnetizationDensity input: "};
  checkCenterCorner(*this, error_tag);

  if (has("grid") && has("dr"))
    throw UniformCommunicateError(error_tag + " cannot define grid and dr.");

  if (has("grid"))
  {
    PosType thisgrid = get<PosType>("grid");
    for (int d = 0; d < DIM; ++d)
    {
      if (thisgrid[d] < 1.0)
        throw UniformCommunicateError(error_tag + " number of grid points must be >=1 in each direction");
    }
  }
  else if (has("dr"))
  {
    //This is the most we can test without knowing the lattice.
    //Correctness determined if dr implies grid with more than 1 point in each direction.
    //Check is performed in calculateDerivedParameters().
    PosType thisdr = get<PosType>("dr");
    for (int d = 0; d < DIM; ++d)
    {
      if (thisdr[d] <= 0.0)
        throw UniformCommunicateError(error_tag + " grid spacing dr must be >= 0 in each direction");
      if (thisdr[d] >= 10.0)
        app_log() << error_tag + " dr larger than 10.0.  Make sure that this grid spacing is intended.\n";
    }
  }
  else
    throw UniformCommunicateError(error_tag + " grid or dr must be specified.");

  if (has("samples"))
  {
    int samps = get<int>("samples");
    if (samps < 1)
      throw UniformCommunicateError(error_tag + " number of samples has to be greater than 0");
  }
}

} //namespace qmcplusplus
