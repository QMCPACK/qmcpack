/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_STRUCTS_STD_H
#define BSPLINE_STRUCTS_STD_H
#include <stdlib.h>

///////////////////////////
// Single precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  Ugrid x_grid;
  BCtype_s xBC;
} UBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
} UBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  size_t coefs_size;
} UBspline_3d_s;


///////////////////////////
// Double precision real //
///////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  size_t coefs_size;
} UBspline_3d_d;



//////////////////////////////
// Single precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  Ugrid x_grid;
  BCtype_c xBC;
} UBspline_1d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_c xBC, yBC;
} UBspline_2d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_float* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_c xBC, yBC, zBC;

} UBspline_3d_c;


//////////////////////////////
// Double precision complex //
//////////////////////////////
typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  Ugrid x_grid;
  BCtype_z xBC;
} UBspline_1d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_z xBC, yBC;
} UBspline_2d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  complex_double* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
} UBspline_3d_z;


#endif
