/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_COPY_H
#define MULTI_BSPLINE_COPY_H

#include "bspline_base.h"
#include "bspline_structs.h"
#include "multi_bspline_structs.h"

#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////              Spline copy functions                 ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/// Create 3D uniform single-precision, real Bspline
  multi_UBspline_3d_s *
  copy_multi_UBspline_3d_s (multi_UBspline_3d_s* spline);

/// Create 3D uniform double-precision, real Bspline
  multi_UBspline_3d_d *
  copy_multi_UBspline_3d_d (multi_UBspline_3d_d* spline);

  /** copy UBspline_3d_d to multi_UBspine_3d_X
   *
   * @param spline target spline
   * @param i state index
   * @param spline_in input spline
   * @param offset offset to match in/out
   * @param N points for each direction
   */
  void copy_UBspline_3d_d(multi_UBspline_3d_d* spline
                          , int i, const UBspline_3d_d* spline_in
                          , const int* offset, const int* N);

  void copy_UBspline_3d_d_s(multi_UBspline_3d_s* spline
                            , int i, const UBspline_3d_d* spline_in
                            , const int* offset, const int* N);

  void copy_UBspline_1d_d(multi_UBspline_1d_d* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N);

  void copy_UBspline_1d_d_s(multi_UBspline_1d_s* multi, int i
      , const UBspline_1d_d* single, const int offset, const int N);

#ifdef __cplusplus
}
#endif

#endif
