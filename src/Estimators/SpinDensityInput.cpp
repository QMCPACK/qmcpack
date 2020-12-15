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

namespace qmcplusplus
{

void SpinDensityInput::readXML(xmlNodePtr cur)
{
  std::string  write_report;
  OhmmsAttributeSet attrib;
  attrib.add(myName_, "name");
  attrib.add(write_report, "report");
  attrib.put(cur);

  bool have_dr     = false;
  bool have_grid   = false;
  bool have_center = false;
  bool have_corner = false;
  bool have_cell   = false;

  PosType dr;
  PosType center;
  Tensor<Real, DIM> axes;

  int test_moves = 0;

  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename((const char*)element->name);
    if (ename == "parameter")
    {
      const XMLAttrString name(element, "name");
      if (name == "dr")
      {
        have_dr = true;
        putContent(dr, element);
      }
      else if (name == "grid")
      {
        have_grid = true;
        putContent(grid_, element);
      }
      else if (name == "corner")
      {
        have_corner = true;
        putContent(corner_, element);
      }
      else if (name == "center")
      {
        have_center = true;
        putContent(center, element);
      }
      else if (name == "cell")
      {
        have_cell = true;
        putContent(axes, element);
      }
      else if (name == "test_moves")
        putContent(test_moves, element);
    }
    element = element->next;
  }

  if (have_dr && have_grid)
  {
    APP_ABORT("SpinDensity::put  dr and grid are provided, this is ambiguous");
  }
  else if (!have_dr && !have_grid)
    APP_ABORT("SpinDensity::put  must provide dr or grid");

  if (have_corner && have_center)
    APP_ABORT("SpinDensity::put  corner and center are provided, this is ambiguous");
  if (have_cell)
  {
    cell_.set(axes);
    if (!have_corner && !have_center)
      APP_ABORT("SpinDensity::put  must provide corner or center");
  }

  if (have_center)
    corner_ = center - cell_.Center;

  if (have_dr)
    for (int d = 0; d < DIM; ++d)
        grid_[d] = (int)std::ceil(std::sqrt(dot(cell_.Rv[d], cell_.Rv[d])) / dr[d]);

  npoints_ = 1;
  for (int d = 0; d < DIM; ++d)
    npoints_ *= grid_[d];
  gdims_[0] = npoints_ / grid_[0];
  for (int d = 1; d < DIM; ++d)
    gdims_[d] = gdims_[d - 1] / grid_[d];
  if (write_report == "yes")
      write_report_ = true;
  else
      write_report_ = false;
  
  // weird legacy stuff
  // if (write_report == "yes")
  //   report("  ");
  // if (test_moves > 0)
  //   test(test_moves, *Ptmp);
}

} // namespace qmcplusplus
