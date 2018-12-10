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
 * @brief Implementation of einspline::Allocator member functions
 *
 * Allocator::Policy is not defined precisely yet but is intended to select 
 * the specialized allocator.
 */
#include "spline2/bspline_allocator.hpp"
#include "spline2/einspline_allocator.h"
#include "einspline/multi_bspline_copy.h"

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

  Allocator::Allocator(): Policy(0){ }

  Allocator::~Allocator()
  {
  }

  multi_UBspline_3d_s* 
    Allocator::allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, int num_splines)
  {
    return einspline_create_multi_UBspline_3d_s(x_grid,y_grid,z_grid,xBC,yBC,zBC,num_splines);
  }

  multi_UBspline_3d_d* 
    Allocator::allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, int num_splines)
  {
    return einspline_create_multi_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,num_splines);
  }

  UBspline_3d_d* Allocator::allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid
      z_grid, BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, double* data)
  {
    return einspline_create_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,data);
  }

  UBspline_3d_s* Allocator::allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid
      z_grid, BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, float* data)
  {
    return einspline_create_UBspline_3d_s(x_grid,y_grid,z_grid,xBC,yBC,zBC,data);
  }

  /** set spline using data in double
   * @param indata raw data
   * @param spline multi spline container 
   * @param i ortbial index
   * 
   * Utility function to use double precision in interpolation and cast to the single-precision spline.
   * This is incomplete and should be integrated with the initialization
   */
  void Allocator::set(double* indata, multi_UBspline_3d_s* spline, int i)
  {
    int BaseOffset[]={0,0,0};
    int BaseN[]={spline->x_grid.num+3, spline->y_grid.num+3, spline->z_grid.num+3};
    BCtype_d xBC_d, yBC_d, zBC_d;
    xBC_d.lCode=spline->xBC.lCode;xBC_d.rCode=spline->xBC.rCode;
    yBC_d.lCode=spline->yBC.lCode;yBC_d.rCode=spline->yBC.rCode;
    zBC_d.lCode=spline->zBC.lCode;zBC_d.rCode=spline->zBC.rCode;
    UBspline_3d_d *singleS=einspline_create_UBspline_3d_d(spline->x_grid,spline->y_grid,spline->z_grid, xBC_d,yBC_d,zBC_d,indata);
    copy(singleS, spline,i,BaseOffset,BaseN);
    einspline_free(singleS->coefs);
    free(singleS);
  }

  void Allocator::set(double* indata, multi_UBspline_3d_d* spline, int i)
  {
    set_multi_UBspline_3d_d(spline, i, indata); 
  }

  void Allocator::set(float* indata, multi_UBspline_3d_s* spline, int i)
  {
    set_multi_UBspline_3d_s(spline, i, indata); 
  }

  void Allocator::copy(multi_UBspline_3d_d* in, multi_UBspline_3d_s* out)
  {
    //Do we ever need this????
    //APP_ABORT("Not done yet with Allocator::copy(in,out)");
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
