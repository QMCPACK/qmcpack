/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_CREATE_H
#define BSPLINE_CREATE_H

#include "bspline_base.h"
#include "bspline_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Spline creation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/////////////////////////////////////
// Uniform, double precision, real //
/////////////////////////////////////
// Create 1D uniform single-precision, real Bspline
  UBspline_1d_d *
  create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data);

// Create 3D uniform single-precision, real Bspline
  UBspline_3d_d *
  create_UBspline_3d_d (Ugrid x_grid,   Ugrid   y_grid,   Ugrid z_grid,
                        BCtype_d  xBC,  BCtype_d   yBC, BCtype_d   zBC,
                        double *data);

  void
  recompute_UBspline_3d_d (UBspline_3d_d* spline, double *data);

#ifdef __cplusplus
}
#endif

#endif
