//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SplineUtils.h"
#include <sstream>
#include "spline2/MultiBsplineBase.hpp"
#include "spline2/MultiBspline1D.hpp"
#include "spline2/einspline_engine.hpp"
#include "spline2/einspline_util.hpp"

namespace qmcplusplus
{
template<typename ST>
bool SplineUtils<ST>::read(MultiBsplineBase<ST>& spline, hdf_archive& h5f)
{
  std::ostringstream o;
  bool success = true;
  int my_index = 0;
  do
  {
    o << "spline_" << my_index;
    einspline_engine<ST, 3> bigtable(spline.getBlock(my_index));
    success = h5f.readEntry(bigtable, o.str());
    my_index++;
  } while (my_index < spline.getNumBlocks() && success);
  return success;
}

template<typename ST>
bool SplineUtils<ST>::write(MultiBsplineBase<ST>& spline, hdf_archive& h5f)
{
  std::ostringstream o;
  bool success = true;
  int my_index = 0;
  do
  {
    o << "spline_" << my_index;
    einspline_engine<ST, 3> bigtable(spline.getBlock(my_index));
    success = h5f.writeEntry(bigtable, o.str());
    my_index++;
  } while (my_index < spline.getNumBlocks() && success);
  return success;
}

template<typename ST>
bool SplineUtils<ST>::read(MultiBspline1D<ST>& spline, hdf_archive& h5f)
{
  einspline_engine<ST, 1> bigtable(*spline.getSplinePtr());
  return h5f.readEntry(bigtable, "radial_spline");
}

template<typename ST>
bool SplineUtils<ST>::write(MultiBspline1D<ST>& spline, hdf_archive& h5f)
{
  einspline_engine<ST, 1> bigtable(*spline.getSplinePtr());
  return h5f.writeEntry(bigtable, "radial_spline");
}

template<typename ST>
void SplineUtils<ST>::gatherv(MultiBsplineBase<ST>& spline,
                              size_t iblock,
                              const std::vector<int>& offset,
                              Communicate& comm)
{
  if (comm.size() == 1)
    return;
  auto& spline_block = spline.getBlock(iblock);
  qmcplusplus::gatherv<ST, 3>(&comm, &spline_block, spline_block.z_stride, offset);
}

template<typename ST>
void SplineUtils<ST>::bcast(MultiBsplineBase<ST>& spline, size_t iblock, Communicate& comm)
{
  if (comm.size() == 1)
    return;
  chunked_bcast(&comm, &spline.getBlock(iblock));
}

template<typename ST>
void SplineUtils<ST>::gatherv(MultiBspline1D<ST>& spline,
                              size_t stride,
                              const std::vector<int>& offset,
                              Communicate& comm)
{
  if (comm.size() == 1)
    return;
  qmcplusplus::gatherv<ST, 1>(&comm, spline.getSplinePtr(), stride, offset);
}

template<typename ST>
void SplineUtils<ST>::bcast(MultiBspline1D<ST>& spline, Communicate& comm)
{
  if (comm.size() == 1)
    return;
  chunked_bcast(&comm, spline.getSplinePtr());
}

template class SplineUtils<float>;
template class SplineUtils<double>;
} // namespace qmcplusplus
