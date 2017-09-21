/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//                                                                         //
//  This program is free software; you can redistribute it and/or modify   //
//  it under the terms of the GNU General Public License as published by   //
//  the Free Software Foundation; either version 2 of the License, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  This program is distributed in the hope that it will be useful,        //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//  GNU General Public License for more details.                           //
//                                                                         //
//  You should have received a copy of the GNU General Public License      //
//  along with this program; if not, write to the Free Software            //
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,                     //
//  Boston, MA  02110-1301  USA                                            //
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
