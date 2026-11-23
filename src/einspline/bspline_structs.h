/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_STRUCTS_STD_H
#define BSPLINE_STRUCTS_STD_H
#include <cstddef>
#include "bspline_base.h"

//////////////////////
// Single precision //
//////////////////////
struct UBspline_1d_s
{
  spline_code spcode;
  type_code tcode;
  float* restrict coefs;
  Ugrid x_grid;
  BCtype_s xBC;
};


struct UBspline_3d_s
{
  spline_code spcode;
  type_code tcode;
  float* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  size_t coefs_size;
};


//////////////////////
// Double precision //
//////////////////////
struct UBspline_1d_d
{
  spline_code spcode;
  type_code tcode;
  double* restrict coefs;
  Ugrid x_grid;
  BCtype_d xBC;
};


struct UBspline_3d_d
{
  spline_code spcode;
  type_code tcode;
  double* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  size_t coefs_size;
};

#endif
