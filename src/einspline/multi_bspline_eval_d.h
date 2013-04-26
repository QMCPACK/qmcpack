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

#ifndef MULTI_BSPLINE_EVAL_D_H
#define MULTI_BSPLINE_EVAL_D_H

#include <math.h>
#include <stdio.h>
#include "multi_bspline_structs.h"

/************************************************************/
/* 1D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_d (const multi_UBspline_1d_d *spline,
                          double x,
                          double* restrict vals);

void
eval_multi_UBspline_1d_d_vg (const multi_UBspline_1d_d *spline,
                             double x,
                             double* restrict vals,
                             double* restrict grads);

void
eval_multi_UBspline_1d_d_vgl (const multi_UBspline_1d_d *spline,
                              double x,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl);

void
eval_multi_UBspline_1d_d_vgh (const multi_UBspline_1d_d *spline,
                              double x,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess);

/************************************************************/
/* 2D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_d (const multi_UBspline_2d_d *spline,
                          double x, double y,
                          double* restrict vals);

void
eval_multi_UBspline_2d_d_vg (const multi_UBspline_2d_d *spline,
                             double x, double y,
                             double* restrict vals,
                             double* restrict grads);

void
eval_multi_UBspline_2d_d_vgl (const multi_UBspline_2d_d *spline,
                              double x, double y,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl);

void
eval_multi_UBspline_2d_d_vgh (const multi_UBspline_2d_d *spline,
                              double x, double y,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess);

/************************************************************/
/* 3D double-precision, real evaulation functions           */
/************************************************************/
void
eval_multi_UBspline_3d_d (const multi_UBspline_3d_d *spline,
                          double x, double y, double z,
                          double* restrict vals);

void
eval_multi_UBspline_3d_d_vg (const multi_UBspline_3d_d *spline,
                             double x, double y, double z,
                             double* restrict vals,
                             double* restrict grads);

void
eval_multi_UBspline_3d_d_vgl (const multi_UBspline_3d_d *spline,
                              double x, double y, double z,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl);

void
eval_multi_UBspline_3d_d_vgh (const multi_UBspline_3d_d *spline,
                              double x, double y, double z,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess);

void
eval_multi_UBspline_3d_d_vghgh (const multi_UBspline_3d_d *spline,
                                double x, double y, double z,
                                double* restrict vals,
                                double* restrict grads,
                                double* restrict hess,
                                double* restrict gradhess);

#endif
