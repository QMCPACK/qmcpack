/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_BASE_H
#define BSPLINE_BASE_H

#define COALLESCED_SIZE 16

#include "config.h"


// Conventions:
// Postfixes:
// s:  single precision real
// d:  double precision real

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Basic type declarations               ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

enum bc_code
{
  PERIODIC,
  DERIV1,
  DERIV2,
  FLAT,
  NATURAL,
  ANTIPERIODIC
};
enum spline_code
{
  U1D,
  U3D,
  MULTI_U1D,
  MULTI_U3D
};
enum type_code
{
  SINGLE_REAL,
  DOUBLE_REAL
};

struct BCtype_s
{
  bc_code lCode, rCode;
  float lVal, rVal;
};

struct BCtype_d
{
  bc_code lCode, rCode;
  double lVal, rVal;
};


struct Ugrid
{
  double start, end;
  int num;

  // private
  double delta, delta_inv;
};

struct Bspline
{
  spline_code sp_code;
  type_code t_code;
  void* restrict coefs;
};

void destroy_Bspline(void* spline);

#endif
