//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file einspline_util.hpp
 * @brief utility functions for i/o and bcast of einspline objects
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_UTILITIES_H
#define QMCPLUSPLUS_EINSPLINE_UTILITIES_H

#include "mpi/mpi_datatype.h"
#include "Message/CommOperators.h"
#include "OhmmsData/FileUtility.h"
#include "hdf/hdf_archive.h"
#include "einspline/multi_bspline_copy.h"
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
{
  chunked_bcast(comm, buffer->coefs, buffer->coefs_size);
}

template<typename ENGT>
inline void gatherv(Communicate* comm, ENGT* buffer, const int ncol, std::vector<int>& offset)
{
  std::vector<int> counts(offset.size() - 1, 0);
  for (size_t ib = 0; ib < counts.size(); ib++)
    counts[ib] = offset[ib + 1] - offset[ib];
  const size_t coef_type_bytes = sizeof(typename bspline_engine_traits<ENGT>::value_type);
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
      comm->gatherv_in_place(buffer->coefs + xs * iz, columntype, counts, offset);
    mpi::free_column_type(columntype);
  }
  else
  {
    const int nrow          = buffer->coefs_size / ncol;
    MPI_Datatype columntype = mpi::construct_column_type(buffer->coefs, nrow, ncol);
    comm->gatherv_in_place(buffer->coefs, columntype, counts, offset);
    mpi::free_column_type(columntype);
  }
}

template<unsigned DIM>
struct dim_traits
{};

// for 3D multi
template<>
struct dim_traits<4>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner->x_grid.num + 3;
    dims[1] = a.spliner->y_grid.num + 3;
    dims[2] = a.spliner->z_grid.num + 3;
    dims[3] = a.spliner->z_stride;
  }
};

// for 1D multi
template<>
struct dim_traits<2>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner->x_grid.num + 2;
    dims[1] = a.spliner->x_stride;
  }
};

/** specialization of h5data_proxy for einspline_engine
   */
template<typename ENGT>
struct h5data_proxy<einspline_engine<ENGT>>
    : public h5_space_type<typename einspline_engine<ENGT>::value_type, einspline_engine<ENGT>::D + 1>
{
  enum
  {
    D = einspline_engine<ENGT>::D
  };
  using value_type = typename einspline_engine<ENGT>::value_type;
  using Base       = h5_space_type<value_type, D + 1>;
  using Base::dims;
  using Base::get_address;
  using data_type = einspline_engine<ENGT>;

  inline h5data_proxy(const data_type& a) { dim_traits<D + 1>::setdim(a, dims); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (ref.spliner)
      return h5d_read(grp, aname, get_address(ref.spliner->coefs), xfer_plist);
    else
      return false;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), Base::rank, dims, get_address(ref.spliner->coefs), xfer_plist);
  }
};

template<typename GT>
void print_grid(GT& grid)
{
  std::cout << grid.start << " " << grid.end << " " << grid.num << " " << grid.delta << " " << grid.delta_inv
            << std::endl;
}

template<typename ENGT>
void print_spliner(ENGT* spline)
{
  std::cout << "xgrid       ";
  print_grid(spline->x_grid);
  std::cout << "ygrid       ";
  print_grid(spline->y_grid);
  std::cout << "zgrid       ";
  print_grid(spline->z_grid);
}


/** engine to create a einspline based on the input einspline
   */
template<typename T>
struct GridConvert
{
  ///number of points in each direction including BC
  int N[3];
  ///offset
  int Offset[3];
  ///number of points of the original grid
  int BaseN[3];
  ///offset of the original grid, always 0
  int BaseOffset[3];
  ///spacing to be added or removed from the buffer
  TinyVector<T, 3> Delta;

  template<typename ENGT1, typename ENGT2, typename PT>
  void create(ENGT1*& out, ENGT2* in, PT& lower, PT& upper, int num)
  {
    using real_type = typename bspline_engine_traits<ENGT1>::real_type;

    Ugrid agrid[3];
    agrid[0] = in->x_grid;
    agrid[1] = in->y_grid;
    agrid[2] = in->z_grid;

    typename bspline_engine_traits<ENGT1>::BCType xyz_bc[3];
    xyz_bc[0].lCode = in->xBC.lCode;
    xyz_bc[0].rCode = in->xBC.rCode;
    xyz_bc[1].lCode = in->yBC.lCode;
    xyz_bc[1].rCode = in->yBC.rCode;
    xyz_bc[2].lCode = in->zBC.lCode;
    xyz_bc[2].rCode = in->zBC.rCode;

    for (int i = 0; i < 3; ++i)
    {
      int ngi        = (int)(lower[i] * agrid[i].delta_inv);
      int ngf        = (int)(upper[i] * agrid[i].delta_inv) + 1;
      agrid[i].num   = std::min(agrid[i].num, ngf - ngi);
      ngf            = agrid[i].num + ngi;
      agrid[i].start = static_cast<real_type>(ngi) * agrid[i].delta;
      agrid[i].end   = static_cast<real_type>(ngf) * agrid[i].delta;

      if (xyz_bc[i].lCode == PERIODIC || xyz_bc[i].lCode == ANTIPERIODIC)
        N[i] = agrid[i].num + 3;
      else
        N[i] = agrid[i].num + 2;

      Delta[i]  = agrid[i].delta;
      Offset[i] = ngi;
    }
    out = einspline::create(out, agrid, xyz_bc, num);
  }
};

namespace einspline
{
template<typename IV>
bool outOfBound(const IV& a)
{
  for (int i = 0; i < a.size(); ++i)
    if (a[i] < 0.0 || a[i] >= 1.0)
      return true;
  return false;
}

template<typename IV>
bool validRange(const IV& low, const IV& up)
{
  bool ok = low[0] < up[0];
  for (int i = 1; i < low.size(); ++i)
    ok &= (low[i] < up[i]);
  return ok;
}
} // namespace einspline

} // namespace qmcplusplus
#endif
