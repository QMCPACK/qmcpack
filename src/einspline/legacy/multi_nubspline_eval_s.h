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

#ifndef MULTI_NBSPLINE_EVAL_S_H
#define MULTI_NBSPLINE_EVAL_S_H

#include "multi_nubspline_structs.h"

/************************************************************/
/* 1D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_1d_s (multi_NUBspline_1d_s *spline,
                          double x,
                          float* restrict vals);

void
eval_multi_NUBspline_1d_s_vg (multi_NUBspline_1d_s *spline,
                             double x,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_NUBspline_1d_s_vgl (multi_NUBspline_1d_s *spline,
                              double x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);


void
eval_multi_NUBspline_1d_s_vgh (multi_NUBspline_1d_s *spline,
                              double x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);

/************************************************************/
/* 2D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_2d_s(multi_NUBspline_2d_s *spline,
                         double x, double y,
                         float* restrict vals);

void
eval_multi_NUBspline_2d_s_vg (multi_NUBspline_2d_s *spline,
                             double x, double y,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_NUBspline_2d_s_vgl (multi_NUBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);

void
eval_multi_NUBspline_2d_s_vgh (multi_NUBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);

/************************************************************/
/* 3D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_3d_s (multi_NUBspline_3d_s *spline,
                          double x, double y, double z,
                          float* restrict vals);

void
eval_multi_NUBspline_3d_s_vg (multi_NUBspline_3d_s *spline,
                             double x, double y, double z,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_NUBspline_3d_s_vgl (multi_NUBspline_3d_s *spline,
                              double x, double y, double z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);

void
eval_multi_NUBspline_3d_s_vgh (multi_NUBspline_3d_s *spline,
                              double x, double y, double z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);
#endif
