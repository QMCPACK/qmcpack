//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "HybridRepCenterOrbitals.h"
#include "hdf/hdf_archive.h"
#include "spline2/SplineUtils.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
template<typename ST>
void AtomicOrbitals<ST>::bcast_tables(Communicate& comm)
{ SplineUtils<ST>::bcast(*SplineInst, comm); }

template<typename ST>
void AtomicOrbitals<ST>::gather_tables(Communicate& comm, const std::vector<int>& offset)
{ SplineUtils<ST>::gatherv(*SplineInst, Npad, offset, comm); }

template<typename ST>
bool AtomicOrbitals<ST>::read_splines(hdf_archive& h5f)
{
  int lmax_in = 0, spline_npoints_in = 0;
  ST spline_radius_in;
  if (!h5f.readEntry(lmax_in, "l_max") || lmax_in != lmax)
    return false;
  if (!h5f.readEntry(spline_radius_in, "spline_radius") || spline_radius_in != spline_radius)
    return false;
  if (!h5f.readEntry(spline_npoints_in, "spline_npoints") || spline_npoints_in != spline_npoints)
    return false;
  return SplineUtils<ST>::read(*SplineInst, h5f);
}

template<typename ST>
bool AtomicOrbitals<ST>::write_splines(hdf_archive& h5f)
{
  bool success = true;
  success      = success && h5f.writeEntry(spline_radius, "spline_radius");
  success      = success && h5f.writeEntry(spline_npoints, "spline_npoints");
  success      = success && h5f.writeEntry(lmax, "l_max");
  success      = success && h5f.writeEntry(center_pos, "position");
  return success && SplineUtils<ST>::write(*SplineInst, h5f);
}

template<typename ST>
void HybridRepCenterOrbitals<ST>::bcast_atomic_tables(Communicate& comm)
{
  for (int ic = 0; ic < AtomicCenters.size(); ic++)
    AtomicCenters[ic].bcast_tables(comm);
}

template<typename ST>
void HybridRepCenterOrbitals<ST>::gather_atomic_tables(Communicate& comm, const std::vector<int>& offset)
{
  if (comm.size() == 1)
    return;
  for (int ic = 0; ic < AtomicCenters.size(); ic++)
    AtomicCenters[ic].gather_tables(comm, offset);
}

template<typename ST>
bool HybridRepCenterOrbitals<ST>::read_atomic_splines(hdf_archive& h5f)
{
  bool success = true;
  size_t ncenter;

  try
  {
    h5f.push("atomic_centers", false);
  }
  catch (...)
  {
    success = false;
  }
  success = success && h5f.readEntry(ncenter, "number_of_centers");
  if (!success)
    return success;
  if (ncenter != AtomicCenters.size())
    success = false;
  // read splines of each center
  for (int ic = 0; ic < AtomicCenters.size(); ic++)
  {
    std::ostringstream gname;
    gname << "center_" << ic;
    try
    {
      h5f.push(gname.str().c_str(), false);
    }
    catch (...)
    {
      success = false;
    }
    success = success && AtomicCenters[ic].read_splines(h5f);
    h5f.pop();
  }
  h5f.pop();
  return success;
}

template<typename ST>
bool HybridRepCenterOrbitals<ST>::write_atomic_splines(hdf_archive& h5f)
{
  bool success = true;
  int ncenter  = AtomicCenters.size();
  try
  {
    h5f.push("atomic_centers", true);
  }
  catch (...)
  {
    success = false;
  }
  success = success && h5f.writeEntry(ncenter, "number_of_centers");
  // write splines of each center
  for (int ic = 0; ic < AtomicCenters.size(); ic++)
  {
    std::ostringstream gname;
    gname << "center_" << ic;
    try
    {
      h5f.push(gname.str().c_str(), true);
    }
    catch (...)
    {
      success = false;
    }
    success = success && AtomicCenters[ic].write_splines(h5f);
    h5f.pop();
  }
  h5f.pop();
  return success;
}

template class AtomicOrbitals<float>;
template class AtomicOrbitals<double>;
template class HybridRepCenterOrbitals<float>;
template class HybridRepCenterOrbitals<double>;
} // namespace qmcplusplus
