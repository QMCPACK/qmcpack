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

#include "OhmmsData/AttributeSet.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

SpinDensityInput::SpinDensityInput(xmlNodePtr node)
{
  readXML(node);
}

void SpinDensityInput::readXML(xmlNodePtr cur)
{
  std::string write_report;
  std::string save_memory;
  OhmmsAttributeSet attrib;
  attrib.add(myName_, "name");
  attrib.add(write_report, "report");
  attrib.add(save_memory, "save_memory");
  attrib.put(cur);

  Tensor<Real, DIM> axes;

  int test_moves = 0;

  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename((const char*)element->name);
    if (ename == "parameter")
    {
      const std::string name(getXMLAttributeValue(element, "name"));
      if (name == "dr")
      {
        have_dr_ = true;
        putContent(dr_, element);
      }
      else if (name == "grid")
      {
        have_grid_ = true;
        putContent(grid_, element);
      }
      else if (name == "corner")
      {
        have_corner_ = true;
        putContent(corner_, element);
      }
      else if (name == "center")
      {
        have_center_ = true;
        putContent(center_, element);
      }
      else if (name == "cell")
      {
        have_cell_ = true;
        putContent(axes, element);
      }
      else if (name == "test_moves")
        putContent(test_moves, element);
    }
    element = element->next;
  }

  if (have_dr_ && have_grid_)
    throw UniformCommunicateError("SpinDensity input dr and grid are provided, this is ambiguous");
  else if (!have_dr_ && !have_grid_)
    throw UniformCommunicateError("SpinDensity input must provide dr or grid");

  if (have_corner_ && have_center_)
    throw UniformCommunicateError("SpinDensity input corner and center are provided, this is ambiguous");

  if (have_cell_)
  {
    cell_.set(axes);
    if (!have_corner_ && !have_center_)
      throw UniformCommunicateError("SpinDensity input must provide corner or center with explicitly defined cell");
  }

  if (write_report == "yes")
    write_report_ = true;
  else
    write_report_ = false;

  if (save_memory == "yes")
    save_memory_ = true;
  else
    save_memory_ = false;

  // weird legacy stuff
  // if (write_report == "yes")
  //   report("  ");
  // if (test_moves > 0)
  //   test(test_moves, *Ptmp);
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
