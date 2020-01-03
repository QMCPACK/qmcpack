//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file bspline_allocator.hpp
 * @brief BsplineAllocator and management classes
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H
#define QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H

#include <simd/Mallocator.hpp>
#include <simd/simd.hpp>
#include <spline2/bspline_traits.hpp>
#include "simd/allocator.hpp"

extern "C"
{
  //forward declaration of the functions
  void find_coefs_1d_d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride);

  void find_coefs_1d_s(Ugrid grid, BCtype_s bc, float* data, intptr_t dstride, float* coefs, intptr_t cstride);
}

namespace qmcplusplus
{
namespace einspline
{
inline void find_coefs_1d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride)
{
  find_coefs_1d_d(grid, bc, data, dstride, coefs, cstride);
}

inline void find_coefs_1d(Ugrid grid, BCtype_s bc, float* data, intptr_t dstride, float* coefs, intptr_t cstride)
{
  find_coefs_1d_s(grid, bc, data, dstride, coefs, cstride);
}
} // namespace einspline

template<typename T, size_t ALIGN = QMC_CLINE, typename ALLOC = Mallocator<T, ALIGN>>
class BsplineAllocator
{
  using SplineType       = typename bspline_traits<T, 3>::SplineType;
  using SingleSplineType = typename bspline_traits<T, 3>::SingleSplineType;
  using BCType           = typename bspline_traits<T, 3>::BCType;
  using real_type        = typename bspline_traits<T, 3>::real_type;

  /// allocator
  ALLOC mAllocator;

public:
  ///default constructor
  BsplineAllocator() = default;
  ///default destructor
  ~BsplineAllocator() = default;
  ///disable copy constructor
  BsplineAllocator(const BsplineAllocator&) = delete;
  ///disable assignement
  BsplineAllocator& operator=(const BsplineAllocator&) = delete;

  template<typename SplineType>
  void destroy(SplineType* spline)
  {
    mAllocator.deallocate(spline->coefs, spline->coefs_size);
    delete (spline);
  }

  ///allocate a multi-bspline structure
  SplineType* allocateMultiBspline(Ugrid x_grid,
                                   Ugrid y_grid,
                                   Ugrid z_grid,
                                   BCType xBC,
                                   BCType yBC,
                                   BCType zBC,
                                   int num_splines);

  ///allocate a UBspline_3d_d, it can be made template to support UBspline_3d_s
  SingleSplineType* allocateUBspline(Ugrid x_grid,
                                     Ugrid y_grid,
                                     Ugrid z_grid,
                                     BCType xBC,
                                     BCType yBC,
                                     BCType zBC,
                                     T* data = nullptr);

  /** copy a UBSpline_3d_X to multi_UBspline_3d_X at i-th band
     * @param single  UBspline_3d_X
     * @param multi target multi_UBspline_3d_X
     * @param i the band index to copy to
     * @param offset starting offset for AoSoA
     * @param N shape of AoSoA
     */
  template<typename UBT, typename MBT>
  void copy(UBT* single, MBT* multi, int i, const int* offset, const int* N);
};

template<typename T, size_t ALIGN, typename ALLOC>
typename BsplineAllocator<T, ALIGN, ALLOC>::SplineType* BsplineAllocator<T, ALIGN, ALLOC>::allocateMultiBspline(
    Ugrid x_grid,
    Ugrid y_grid,
    Ugrid z_grid,
    BCType xBC,
    BCType yBC,
    BCType zBC,
    int num_splines)
{
  // Create new spline
  SplineType* spline = new SplineType;

  spline->spcode      = bspline_traits<T, 3>::spcode;
  spline->tcode       = bspline_traits<T, 3>::tcode;
  spline->xBC         = xBC;
  spline->yBC         = yBC;
  spline->zBC         = zBC;
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx + 3;
  else
    Nx = Mx + 2;
  x_grid.delta     = (x_grid.end - x_grid.start) / (double)(Nx - 3);
  x_grid.delta_inv = 1.0 / x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My + 3;
  else
    Ny = My + 2;
  y_grid.delta     = (y_grid.end - y_grid.start) / (double)(Ny - 3);
  y_grid.delta_inv = 1.0 / y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
    Nz = Mz + 3;
  else
    Nz = Mz + 2;
  z_grid.delta     = (z_grid.end - z_grid.start) / (double)(Nz - 3);
  z_grid.delta_inv = 1.0 / z_grid.delta;
  spline->z_grid   = z_grid;

  const int N = getAlignedSize<real_type, ALIGN>(num_splines);

  spline->x_stride = (size_t)Ny * (size_t)Nz * (size_t)N;
  spline->y_stride = Nz * N;
  spline->z_stride = N;

  spline->coefs_size = (size_t)Nx * spline->x_stride;
  spline->coefs      = mAllocator.allocate(spline->coefs_size);

  return spline;
}

template<typename T, size_t ALIGN, typename ALLOC>
typename BsplineAllocator<T, ALIGN, ALLOC>::SingleSplineType* BsplineAllocator<T, ALIGN, ALLOC>::allocateUBspline(
    Ugrid x_grid,
    Ugrid y_grid,
    Ugrid z_grid,
    BCType xBC,
    BCType yBC,
    BCType zBC,
    T* data)
{
  // Create new spline
  SingleSplineType* spline = new SingleSplineType;
  spline->spcode           = bspline_traits<T, 3>::single_spcode;
  spline->tcode            = bspline_traits<T, 3>::tcode;
  spline->xBC              = xBC;
  spline->yBC              = yBC;
  spline->zBC              = zBC;

  // Setup internal variables
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx + 3;
  else
    Nx = Mx + 2;
  x_grid.delta     = (x_grid.end - x_grid.start) / (double)(Nx - 3);
  x_grid.delta_inv = 1.0 / x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My + 3;
  else
    Ny = My + 2;
  y_grid.delta     = (y_grid.end - y_grid.start) / (double)(Ny - 3);
  y_grid.delta_inv = 1.0 / y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
    Nz = Mz + 3;
  else
    Nz = Mz + 2;
  z_grid.delta     = (z_grid.end - z_grid.start) / (double)(Nz - 3);
  z_grid.delta_inv = 1.0 / z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny * Nz;
  spline->y_stride = Nz;

  spline->coefs_size = (size_t)Nx * (size_t)Ny * (size_t)Nz;

  spline->coefs = mAllocator.allocate(spline->coefs_size);

  if (data != NULL) // only data is provided
  {
// First, solve in the X-direction
#pragma omp parallel for
    for (int iy = 0; iy < My; iy++)
      for (int iz = 0; iz < Mz; iz++)
      {
        intptr_t doffset = iy * Mz + iz;
        intptr_t coffset = iy * Nz + iz;
        einspline::find_coefs_1d(spline->x_grid, xBC, data + doffset, My * Mz, spline->coefs + coffset, Ny * Nz);
      }

// Now, solve in the Y-direction
#pragma omp parallel for
    for (intptr_t ix = 0; ix < Nx; ix++)
      for (intptr_t iz = 0; iz < Nz; iz++)
      {
        intptr_t doffset = ix * Ny * Nz + iz;
        intptr_t coffset = ix * Ny * Nz + iz;
        einspline::find_coefs_1d(spline->y_grid, yBC, spline->coefs + doffset, Nz, spline->coefs + coffset, Nz);
      }

// Now, solve in the Z-direction
#pragma omp parallel for
    for (intptr_t ix = 0; ix < Nx; ix++)
      for (intptr_t iy = 0; iy < Ny; iy++)
      {
        intptr_t doffset = (ix * Ny + iy) * Nz;
        intptr_t coffset = (ix * Ny + iy) * Nz;
        einspline::find_coefs_1d(spline->z_grid, zBC, spline->coefs + doffset, 1, spline->coefs + coffset, 1);
      }
  }
  return spline;
}

template<typename T, size_t ALIGN, typename ALLOC>
template<typename UBT, typename MBT>
void BsplineAllocator<T, ALIGN, ALLOC>::copy(UBT* single, MBT* multi, int i, const int* offset, const int* N)
{
  typedef typename bspline_type<MBT>::value_type out_type;
  typedef typename bspline_type<UBT>::value_type in_type;
  intptr_t x_stride_in  = single->x_stride;
  intptr_t y_stride_in  = single->y_stride;
  intptr_t x_stride_out = multi->x_stride;
  intptr_t y_stride_out = multi->y_stride;
  intptr_t z_stride_out = multi->z_stride;
  intptr_t offset0      = static_cast<intptr_t>(offset[0]);
  intptr_t offset1      = static_cast<intptr_t>(offset[1]);
  intptr_t offset2      = static_cast<intptr_t>(offset[2]);
  const intptr_t istart = static_cast<intptr_t>(i);
  const intptr_t n0 = N[0], n1 = N[1], n2 = N[2];
  for (intptr_t ix = 0; ix < n0; ++ix)
    for (intptr_t iy = 0; iy < n1; ++iy)
    {
      out_type* restrict out = multi->coefs + ix * x_stride_out + iy * y_stride_out + istart;
      const in_type* restrict in =
          single->coefs + (ix + offset0) * x_stride_in + (iy + offset1) * y_stride_in + offset2;
      for (intptr_t iz = 0; iz < n2; ++iz)
      {
        out[iz * z_stride_out] = static_cast<out_type>(in[iz]);
      }
    }
}

} // namespace qmcplusplus
#endif
