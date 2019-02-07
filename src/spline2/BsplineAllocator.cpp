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
/** @file einspline_allocator.cpp
 * @brief Implementation of einspline::BsplineAllocator member functions
 *
 * BsplineAllocator::Policy is not defined precisely yet but is intended to select
 * the specialized allocator.
 */
#include "spline2/BsplineAllocator.hpp"
#include "spline2/einspline_allocator.h"
#include "einspline/multi_bspline_copy.h"
#include "einspline/multi_bspline_create.h"

extern "C" {

  void set_multi_UBspline_3d_d (multi_UBspline_3d_d *spline, int spline_num, double *data);

  void set_multi_UBspline_3d_s (multi_UBspline_3d_s *spline, int spline_num, float *data);

  multi_UBspline_3d_s*
    einspline_create_multi_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, int num_splines);

  UBspline_3d_s*
    einspline_create_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, float *data);

  multi_UBspline_3d_d*
    einspline_create_multi_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, int num_splines);

  UBspline_3d_d*
    einspline_create_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, double *data);

}

namespace qmcplusplus { namespace einspline {

  BsplineAllocator::BsplineAllocator(): Policy(0){ }

  BsplineAllocator::~BsplineAllocator()
  {
  }

  multi_UBspline_3d_s*
    BsplineAllocator::allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, int num_splines)
  {
    return einspline_create_multi_UBspline_3d_s(x_grid,y_grid,z_grid,xBC,yBC,zBC,num_splines);
  }

  multi_UBspline_3d_d*
    BsplineAllocator::allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, int num_splines)
  {
    return einspline_create_multi_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,num_splines);
  }

  UBspline_3d_d* BsplineAllocator::allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid
      z_grid, BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, double* data)
  {
    return einspline_create_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,data);
  }

  UBspline_3d_s* BsplineAllocator::allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid
      z_grid, BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, float* data)
  {
    return einspline_create_UBspline_3d_s(x_grid,y_grid,z_grid,xBC,yBC,zBC,data);
  }

  // 1D spline interface to einspline routines.
  void set(multi_UBspline_1d_d* spline, int i, UBspline_1d_d* spline_in,
       const int offset, const int N)
  {
    copy_UBspline_1d_d(spline, i, spline_in, offset, N);
  }

  void set(multi_UBspline_1d_s* spline, int i, UBspline_1d_d* spline_in,
       const int offset, const int N)
  {
    copy_UBspline_1d_d_s(spline, i, spline_in, offset, N);
  }

  multi_UBspline_1d_d* create(multi_UBspline_1d_d* s, Ugrid& grid, BCtype_d& bc, int num_splines)
  {
    multi_UBspline_1d_d* newspline=create_multi_UBspline_1d_d(grid, bc, num_splines);
    std::fill(newspline->coefs, newspline->coefs+newspline->coefs_size, 0.0);
    return newspline;
  }

  multi_UBspline_1d_s* create(multi_UBspline_1d_s* s, Ugrid& grid, BCtype_s& bc, int num_splines)
  {
    multi_UBspline_1d_s* newspline=create_multi_UBspline_1d_s(grid, bc, num_splines);
    std::fill(newspline->coefs, newspline->coefs+newspline->coefs_size, 0.0f);
    return newspline;
  }

} }
