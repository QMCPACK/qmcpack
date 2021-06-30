/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_NUBSPLINE_STRUCTS_STD_H
#define MULTI_NUBSPLINE_STRUCTS_STD_H

#include <stdint.h>
#include "bspline_base.h"
#include "nubasis.h"

///////////////////////////
// Single precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride;
  BCtype_s xBC;
  int num_splines;
  NUgrid  *restrict x_grid;
  NUBasis *restrict x_basis;
} multi_NUBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride;
  BCtype_s xBC, yBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid;
  NUBasis *restrict x_basis, *restrict y_basis;
} multi_NUBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  BCtype_s xBC, yBC, zBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid,  *restrict z_grid;
  NUBasis *restrict x_basis, *restrict y_basis, *restrict z_basis;
} multi_NUBspline_3d_s;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride;
  BCtype_d xBC;
  int num_splines;
  NUgrid  *restrict x_grid;
  NUBasis *restrict x_basis;
} multi_NUBspline_1d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride, y_stride;
  BCtype_d xBC, yBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid;
  NUBasis *restrict x_basis, *restrict y_basis;
} multi_NUBspline_2d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  BCtype_d xBC, yBC, zBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid,  *restrict z_grid;
  NUBasis *restrict x_basis, *restrict y_basis, *restrict z_basis;
} multi_NUBspline_3d_d;



//////////////////////////////
// Single precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride;
  BCtype_c xBC;
  int num_splines;
  NUgrid  *restrict x_grid;
  NUBasis *restrict x_basis;
} multi_NUBspline_1d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride, y_stride;
  BCtype_c xBC, yBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid;
  NUBasis *restrict x_basis, *restrict y_basis;
} multi_NUBspline_2d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  BCtype_c xBC, yBC, zBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid,  *restrict z_grid;
  NUBasis *restrict x_basis, *restrict y_basis, *restrict z_basis;
} multi_NUBspline_3d_c;


//////////////////////////////
// Double precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride;
  BCtype_z xBC;
  int num_splines;
  NUgrid  *restrict x_grid;
  NUBasis *restrict x_basis;
} multi_NUBspline_1d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride, y_stride;
  BCtype_z xBC, yBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid;
  NUBasis *restrict x_basis, *restrict y_basis;
} multi_NUBspline_2d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  BCtype_z xBC, yBC, zBC;
  int num_splines;
  NUgrid  *restrict x_grid,  *restrict y_grid,  *restrict z_grid;
  NUBasis *restrict x_basis, *restrict y_basis, *restrict z_basis;
} multi_NUBspline_3d_z;


#endif
