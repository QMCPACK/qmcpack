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

typedef enum
{
  PERIODIC,
  DERIV1,
  DERIV2,
  FLAT,
  NATURAL,
  ANTIPERIODIC
} bc_code;
typedef enum
{
  U1D,
  U3D,
  MULTI_U1D,
  MULTI_U3D
} spline_code;
typedef enum
{
  SINGLE_REAL,
  DOUBLE_REAL
} type_code;

typedef struct
{
  bc_code lCode, rCode;
  float lVal, rVal;
} BCtype_s;

typedef struct
{
  bc_code lCode, rCode;
  double lVal, rVal;
} BCtype_d;


typedef struct
{
  double start, end;
  int num;

  // private
  double delta, delta_inv;
} Ugrid;

typedef struct
{
  spline_code sp_code;
  type_code t_code;
  void* restrict coefs;
} Bspline;

#ifdef __cplusplus
extern "C"
#endif
    void destroy_Bspline(void* spline);

#endif
