/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007-2010 Kenneth P. Esler, Jr.                          //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_BASE_CUDA_H
#define BSPLINE_BASE_CUDA_H

#include <cuda.h>

#if defined(CUDA_VERSION) && (CUDA_VERSION < 3000) /* 3.0 */
typedef struct
{
  double x,y,z;
} double3;

typedef struct
{
  double x,y,z,w;
} double4;
#endif

#endif
