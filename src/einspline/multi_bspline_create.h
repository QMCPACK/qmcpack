/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_CREATE_H
#define MULTI_BSPLINE_CREATE_H

#include "bspline_base.h"
#include "multi_bspline_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Spline creation functions             ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

///////////////////////////////
// Uniform, single precision //
///////////////////////////////

// Create 1D uniform single-precision Bspline
  multi_UBspline_1d_s *
  create_multi_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, int num_splines);

///////////////////////////////
// Uniform, double precision //
///////////////////////////////

// Create 1D uniform single-precision Bspline
  multi_UBspline_1d_d *
  create_multi_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, int num_splines);

#ifdef __cplusplus
}
#endif

#endif
