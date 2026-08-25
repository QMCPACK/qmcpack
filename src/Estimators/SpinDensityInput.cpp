//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SpinDensity.cpp
//////////////////////////////////////////////////////////////////////////////////////
#include "SpinDensityInput.h"

#include <cmath>

#include "EstimatorInput.h"

namespace qmcplusplus
{

SpinDensityInput::SpinDensityInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;

  setIfInInput(name_, "name");
  setIfInInput(write_report_, "report");
  setIfInInput(save_memory_, "save_memory");
  have_dr_     = setIfInInput(dr_, "dr");
  have_corner_ = setIfInInput(corner_, "corner");
  have_center_ = setIfInInput(center_, "center");

  PosType input_grid;
  have_grid_ = setIfInInput(input_grid, "grid");
  if (have_grid_)
    for (int d = 0; d < DIM; ++d)
      grid_[d] = static_cast<int>(input_grid[d]);

  have_cell_ = input_section_.has("cell");
  if (have_cell_)
  {
    const std::vector<Real> cell_values = input_section_.get<std::vector<Real>>("cell");
    Tensor<Real, DIM> axes;
    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        axes(i, j) = cell_values[i * DIM + j];
    cell_.set(axes);
  }
}

void SpinDensityInput::SpinDensityInputSection::checkParticularValidity()
{
  using namespace estimatorinput;
  const std::string error_tag{"SpinDensity input: "};

  checkCenterCorner(*this, error_tag);

  if (has("dr") && has("grid"))
    throw UniformCommunicateError(error_tag + "dr and grid are provided, this is ambiguous");
  if (!has("dr") && !has("grid"))
    throw UniformCommunicateError(error_tag + "must provide dr or grid");

  if (has("cell"))
  {
    if (get<std::vector<Real>>("cell").size() != DIM * DIM)
      throw UniformCommunicateError(error_tag + "cell must contain a " + std::to_string(DIM) + " by " +
                                    std::to_string(DIM) + " matrix");
    if (!has("corner") && !has("center"))
      throw UniformCommunicateError(error_tag + "must provide corner or center with explicitly defined cell");
  }
}

SpinDensityInput::DerivedParameters SpinDensityInput::calculateDerivedParameters(const Lattice& lattice) const
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

} // namespace qmcplusplus
