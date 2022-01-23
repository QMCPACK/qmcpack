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

#include "spline2/bspline_traits.hpp"
#include "CPU/SIMD/aligned_allocator.hpp"

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

template<typename T,
         typename COEFS_ALLOC         = aligned_allocator<T>,
         typename MULTI_SPLINE_ALLOC  = aligned_allocator<typename bspline_traits<T, 3>::SplineType>,
         typename SINGLE_SPLINE_ALLOC = aligned_allocator<typename bspline_traits<T, 3>::SingleSplineType>>
class BsplineAllocator
{
  using SplineType       = typename bspline_traits<T, 3>::SplineType;
  using SingleSplineType = typename bspline_traits<T, 3>::SingleSplineType;
  using BCType           = typename bspline_traits<T, 3>::BCType;
  using real_type        = typename bspline_traits<T, 3>::real_type;

  /// allocators
  COEFS_ALLOC coefs_allocator;
  MULTI_SPLINE_ALLOC multi_spline_allocator;
  SINGLE_SPLINE_ALLOC single_spline_allocator;

public:
  ///default constructor
  BsplineAllocator() = default;
  ///default destructor
  ~BsplineAllocator() = default;
  ///disable copy constructor
  BsplineAllocator(const BsplineAllocator&) = delete;
  ///disable assignement
  BsplineAllocator& operator=(const BsplineAllocator&) = delete;

  void destroy(SplineType* spline)
  {
    coefs_allocator.deallocate(spline->coefs, spline->coefs_size);
    multi_spline_allocator.deallocate(spline, 1);
  }

  void destroy(SingleSplineType* spline)
  {
    coefs_allocator.deallocate(spline->coefs, spline->coefs_size);
    single_spline_allocator.deallocate(spline, 1);
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

template<typename T, typename COEFS_ALLOC, typename MULTI_SPLINE_ALLOC, typename SINGLE_SPLINE_ALLOC>
typename BsplineAllocator<T, COEFS_ALLOC, MULTI_SPLINE_ALLOC, SINGLE_SPLINE_ALLOC>::SplineType* BsplineAllocator<
    T,
    COEFS_ALLOC,
    MULTI_SPLINE_ALLOC,
    SINGLE_SPLINE_ALLOC>::allocateMultiBspline(Ugrid x_grid,
                                               Ugrid y_grid,
                                               Ugrid z_grid,
                                               BCType xBC,
                                               BCType yBC,
                                               BCType zBC,
                                               int num_splines)
{
  // Create new spline
  SplineType* spline = multi_spline_allocator.allocate(1);

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

  const int N = getAlignedSize<real_type, COEFS_ALLOC::alignment>(num_splines);

  spline->x_stride = (size_t)Ny * (size_t)Nz * (size_t)N;
  spline->y_stride = Nz * N;
  spline->z_stride = N;

  spline->coefs_size = (size_t)Nx * spline->x_stride;
  spline->coefs      = coefs_allocator.allocate(spline->coefs_size);

  return spline;
}

template<typename T, typename COEFS_ALLOC, typename MULTI_SPLINE_ALLOC, typename SINGLE_SPLINE_ALLOC>
typename BsplineAllocator<T, COEFS_ALLOC, MULTI_SPLINE_ALLOC, SINGLE_SPLINE_ALLOC>::SingleSplineType* BsplineAllocator<
    T,
    COEFS_ALLOC,
    MULTI_SPLINE_ALLOC,
    SINGLE_SPLINE_ALLOC>::allocateUBspline(Ugrid x_grid,
                                           Ugrid y_grid,
                                           Ugrid z_grid,
                                           BCType xBC,
                                           BCType yBC,
                                           BCType zBC,
                                           T* data)
{
  // Create new spline
  SingleSplineType* spline = single_spline_allocator.allocate(1);
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

  spline->coefs = coefs_allocator.allocate(spline->coefs_size);

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

template<typename T, typename COEFS_ALLOC, typename MULTI_SPLINE_ALLOC, typename SINGLE_SPLINE_ALLOC>
template<typename UBT, typename MBT>
void BsplineAllocator<T, COEFS_ALLOC, MULTI_SPLINE_ALLOC, SINGLE_SPLINE_ALLOC>::copy(UBT* single,
                                                                                     MBT* multi,
                                                                                     int i,
                                                                                     const int* offset,
                                                                                     const int* N)
{
  using out_type        = typename bspline_type<MBT>::value_type;
  using in_type         = typename bspline_type<UBT>::value_type;
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
