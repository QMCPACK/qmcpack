//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file MultiBsplineAllocator.hpp
 * @brief MultiBsplineAllocator classes
 */
#ifndef QMCPLUSPLUS_MULTI_BSPLINE_ALLOCATOR_H
#define QMCPLUSPLUS_MULTI_BSPLINE_ALLOCATOR_H

#include "spline2/bspline_traits.hpp"
#include "CPU/SIMD/aligned_allocator.hpp"

namespace qmcplusplus
{
template<typename T, typename ALLOC = aligned_allocator<T>>
class MultiBsplineAllocator
{
  using SplineType       = typename bspline_traits<T, 3>::SplineType;
  using SingleSplineType = typename bspline_traits<T, 3>::SingleSplineType;
  using BCType           = typename bspline_traits<T, 3>::BCType;
  using real_type        = typename bspline_traits<T, 3>::real_type;

  /// allocators
  ALLOC coefs_allocator;
  using MultiAlloc = typename std::allocator_traits<ALLOC>::template rebind_alloc<SplineType>;
  MultiAlloc multi_spline_allocator;

public:
  ///default constructor
  MultiBsplineAllocator() = default;
  ///default destructor
  ~MultiBsplineAllocator() = default;
  ///disable copy constructor
  MultiBsplineAllocator(const MultiBsplineAllocator&) = delete;
  ///disable assignement
  MultiBsplineAllocator& operator=(const MultiBsplineAllocator&) = delete;

  void destroy(SplineType* spline)
  {
    coefs_allocator.deallocate(spline->coefs, spline->coefs_size);
    multi_spline_allocator.deallocate(spline, 1);
  }

  ///allocate a multi-bspline structure
  SplineType* allocateMultiBspline(Ugrid x_grid,
                                   Ugrid y_grid,
                                   Ugrid z_grid,
                                   const BCType& xBC,
                                   const BCType& yBC,
                                   const BCType& zBC,
                                   int num_splines);
};

template<typename T, typename ALLOC>
typename MultiBsplineAllocator<T, ALLOC>::SplineType* MultiBsplineAllocator<T, ALLOC>::allocateMultiBspline(
    Ugrid x_grid,
    Ugrid y_grid,
    Ugrid z_grid,
    const BCType& xBC,
    const BCType& yBC,
    const BCType& zBC,
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

  const int N = getAlignedSize<real_type, ALLOC::alignment>(num_splines);

  spline->x_stride = (size_t)Ny * (size_t)Nz * (size_t)N;
  spline->y_stride = Nz * N;
  spline->z_stride = N;

  spline->coefs_size = (size_t)Nx * spline->x_stride;
  spline->coefs      = coefs_allocator.allocate(spline->coefs_size);

  return spline;
}

} // namespace qmcplusplus
#endif
