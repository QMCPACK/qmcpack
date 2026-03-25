//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file einspline_util.hpp
 * @brief utility functions for bcast gatherv of einspline objects
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_UTILITIES_H
#define QMCPLUSPLUS_EINSPLINE_UTILITIES_H

#include "Message/CommOperators.h"
#include "mpi/mpi_datatype.h"
#include "Host/OutputManager.h"
#include "bspline_traits.hpp"
#include <limits>

namespace qmcplusplus
{
///handles i/o and bcast, testing for now
template<typename T>
inline void chunked_bcast(Communicate* comm, T* buffer, size_t ntot)
{
  if (comm->size() == 1)
    return;

  size_t chunk_size = (1 << 30) / sizeof(T); //256 MB
  int n             = static_cast<int>(ntot / chunk_size);

  size_t offset = 0;
  for (int i = 0; i < n; ++i, offset += chunk_size)
  {
    comm->bcast(buffer + offset, static_cast<int>(chunk_size));
  }

  if (offset < ntot)
  {
    comm->bcast(buffer + offset, static_cast<int>(ntot - offset));
  }
}

template<typename ENGT>
inline void chunked_bcast(Communicate* comm, ENGT* buffer)
{ chunked_bcast(comm, buffer->coefs, buffer->coefs_size); }

template<typename T, unsigned D>
inline void gatherv(Communicate* comm,
                    typename bspline_traits<T, D>::SplineType* buffer,
                    const int ncol,
                    const std::vector<int>& offset)
{
  std::vector<int> counts(offset.size() - 1, 0);
  for (size_t ib = 0; ib < counts.size(); ib++)
    counts[ib] = offset[ib + 1] - offset[ib];
  const auto& counts_const(counts);
  const size_t coef_type_bytes = sizeof(T);
  if (buffer->coefs_size * coef_type_bytes > std::numeric_limits<int>::max())
  {
    // Some MPI libraries have problems when message sizes exceed range of integer (2^31-1)
    // Perform the gatherv in columns to reduce risk
    const size_t xs = buffer->x_stride;
    if (xs * coef_type_bytes >= std::numeric_limits<int>::max())
      app_warning() << "Large single message even after splitting by the number of grid points in x direction! "
                    << "Some MPI libraries may not work!" << std::endl;
    const size_t nx         = buffer->coefs_size / xs;
    const int nrow          = buffer->coefs_size / (ncol * nx);
    MPI_Datatype columntype = mpi::construct_column_type(buffer->coefs, nrow, ncol);
    for (size_t iz = 0; iz < nx; iz++)
      comm->gatherv_in_place(buffer->coefs + xs * iz, columntype, counts_const, offset);
    mpi::free_column_type(columntype);
  }
  else
  {
    const int nrow          = buffer->coefs_size / ncol;
    MPI_Datatype columntype = mpi::construct_column_type(buffer->coefs, nrow, ncol);
    comm->gatherv_in_place(buffer->coefs, columntype, counts_const, offset);
    mpi::free_column_type(columntype);
  }
}

} // namespace qmcplusplus
#endif
