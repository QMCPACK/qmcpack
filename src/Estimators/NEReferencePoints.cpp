//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak. doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCHamiltonian/ReferencePoints.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "NEReferencePoints.h"
#include "Utilities/string_utils.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCHamiltonians/ObservableHelper.h"

namespace qmcplusplus
{

NEReferencePoints::NEReferencePoints(const ReferencePointsInput& rp_input,
                                     const ParticleSet& pset,
                                     RefVector<ParticleSet>& ref_psets)
    : input_(rp_input)
{
  processParticleSets(pset, ref_psets);
  for (int i = 0; i < OHMMS_DIM; i++)
    for (int d = 0; d < OHMMS_DIM; d++)
      axes(d, i) = pset.getLattice().a(i)[d];
  Axes crd;
  // no need to handle error here rp_input will have a valid value for coord_form
  switch (input_.get_coord_form())
  {
  case Coord::CELL:
    crd = axes;
    break;
  case Coord::CARTESIAN:
    for (int i = 0; i < OHMMS_DIM; i++)
      for (int d = 0; d < OHMMS_DIM; d++)
        if (d == i)
          crd(i, i) = 1.0;
        else
          crd(d, i) = 0.0;
    break;
  }

  for (const auto& [key, value] : input_.get_points())
    points_[key] = dot(crd, value);
}

void NEReferencePoints::processParticleSets(const ParticleSet& P, RefVector<ParticleSet>& Psets)
{
  //get axes and origin information from the ParticleSet
  points_["zero"] = 0 * P.getLattice().a(0);
  points_["a1"]   = P.getLattice().a(0);
  points_["a2"]   = P.getLattice().a(1);
  points_["a3"]   = P.getLattice().a(2);
  //points_["center"]= .5*(P.getLattice().a(0)+P.getLattice().a(1)+P.Lattice.a(2))
  //set points_ on face centers
  points_["f1p"] = points_["zero"] + .5 * points_["a1"];
  points_["f1m"] = points_["zero"] - .5 * points_["a1"];
  points_["f2p"] = points_["zero"] + .5 * points_["a2"];
  points_["f2m"] = points_["zero"] - .5 * points_["a2"];
  points_["f3p"] = points_["zero"] + .5 * points_["a3"];
  points_["f3m"] = points_["zero"] - .5 * points_["a3"];
  //set points_ on cell corners
  points_["cmmm"] = points_["zero"] + .5 * (-1 * points_["a1"] - points_["a2"] - points_["a3"]);
  points_["cpmm"] = points_["zero"] + .5 * (points_["a1"] - points_["a2"] - points_["a3"]);
  points_["cmpm"] = points_["zero"] + .5 * (-1 * points_["a1"] + points_["a2"] - points_["a3"]);
  points_["cmmp"] = points_["zero"] + .5 * (-1 * points_["a1"] - points_["a2"] + points_["a3"]);
  points_["cmpp"] = points_["zero"] + .5 * (-1 * points_["a1"] + points_["a2"] + points_["a3"]);
  points_["cpmp"] = points_["zero"] + .5 * (points_["a1"] - points_["a2"] + points_["a3"]);
  points_["cppm"] = points_["zero"] + .5 * (points_["a1"] + points_["a2"] - points_["a3"]);
  points_["cppp"] = points_["zero"] + .5 * (points_["a1"] + points_["a2"] + points_["a3"]);
  //get points from requested particle sets
  int cshift = 1;
  for (ParticleSet& pset : Psets)
  {
    for (int p = 0; p < pset.getTotalNum(); p++)
    {
      std::stringstream ss;
      ss << p + cshift;
      points_[pset.getName() + ss.str()] = pset.R[p];
    }
  }
}

void NEReferencePoints::write_description(std::ostream& os, const std::string& indent) const
{
  os << indent + "reference_points" << std::endl;
  std::map<std::string, Point>::const_iterator it, end = points_.end();
  for (it = points_.begin(); it != end; ++it)
  {
    os << indent + "  " << it->first << ": " << it->second << std::endl;
  }
  os << indent + "end reference_points" << std::endl;
  return;
}

void NEReferencePoints::write(hdf_archive& file) const
{
  file.push(std::string_view("reference_points"));
  for (auto it = points_.cbegin(); it != points_.cend(); ++it)
    file.write(const_cast<Point&>(it->second), it->first);
  file.pop();
}

std::ostream& operator<<(std::ostream& out, const NEReferencePoints& rhs)
{
  rhs.write_description(out, "");
  return out;
}

} // namespace qmcplusplus
