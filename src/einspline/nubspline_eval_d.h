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

#ifndef NUBSPLINE_EVAL_D_H
#define NUBSPLINE_EVAL_D_H

#include <math.h>
#include <stdio.h>
#include "nubspline_structs.h"

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

void
eval_NUBspline_1d_d (NUBspline_1d_d * restrict spline,
                     double x, double* restrict val);

void
eval_NUBspline_1d_d_vg (NUBspline_1d_d * restrict spline, double x,
                        double* restrict val, double* restrict grad);

void
eval_NUBspline_1d_d_vgl (NUBspline_1d_d * restrict spline, double x,
                         double* restrict val, double* restrict grad,
                         double* restrict lapl);

void
eval_NUBspline_1d_d_vgh (NUBspline_1d_d * restrict spline, double x,
                         double* restrict val, double* restrict grad,
                         double* restrict hess);

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

void
eval_NUBspline_2d_d (NUBspline_2d_d * restrict spline,
                     double x, double y, double* restrict val);

void
eval_NUBspline_2d_d_vg (NUBspline_2d_d * restrict spline,
                        double x, double y,
                        double* restrict val, double* restrict grad);

void
eval_NUBspline_2d_d_vgl (NUBspline_2d_d * restrict spline,
                         double x, double y, double* restrict val,
                         double* restrict grad, double* restrict lapl);

void
eval_NUBspline_2d_d_vgh (NUBspline_2d_d * restrict spline,
                         double x, double y, double* restrict val,
                         double* restrict grad, double* restrict hess);

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

void
eval_NUBspline_3d_d (NUBspline_3d_d * restrict spline,
                     double x, double y, double z,
                     double* restrict val);

void
eval_NUBspline_3d_d_vg (NUBspline_3d_d * restrict spline,
                        double x, double y, double z,
                        double* restrict val, double* restrict grad);

void
eval_NUBspline_3d_d_vgl (NUBspline_3d_d * restrict spline,
                         double x, double y, double z,
                         double* restrict val, double* restrict grad, double* restrict lapl);

void
eval_NUBspline_3d_d_vgh (NUBspline_3d_d * restrict spline,
                         double x, double y, double z,
                         double* restrict val, double* restrict grad, double* restrict hess);

#endif
