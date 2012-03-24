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

#ifndef MULTI_NUBSPLINE_EVAL_C_H
#define MULTI_NUBSPLINE_EVAL_C_H

#include <math.h>
#include <stdio.h>
#include "multi_nubspline_structs.h"

/************************************************************/
/* 1D float-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_1d_c (multi_NUBspline_1d_c *spline,
			   double x,
			   complex_float* restrict vals);

void
eval_multi_NUBspline_1d_c_vg (multi_NUBspline_1d_c *spline,
			      double x,
			      complex_float* restrict vals,
			      complex_float* restrict grads);


void
eval_multi_NUBspline_1d_c_vgl (multi_NUBspline_1d_c *spline,
			       double x,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict lapl);


void
eval_multi_NUBspline_1d_c_vgh (multi_NUBspline_1d_c *spline,
			       double x,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict hess);


/************************************************************/
/* 2D float-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_2d_c (multi_NUBspline_2d_c *spline,
			   double x, double y,
			   complex_float* restrict vals);

void
eval_multi_NUBspline_2d_c_vg (multi_NUBspline_2d_c *spline,
			      double x, double y,
			      complex_float* restrict vals,
			      complex_float* restrict grads);

void
eval_multi_NUBspline_2d_c_vgl (multi_NUBspline_2d_c *spline,
			       double x, double y,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict lapl);

void
eval_multi_NUBspline_2d_c_vgh (multi_NUBspline_2d_c *spline,
			       double x, double y,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict hess);


/************************************************************/
/* 3D float-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_3d_c (multi_NUBspline_3d_c *spline,
			   double x, double y, double z,
			   complex_float* restrict vals);

void
eval_multi_NUBspline_3d_c_vg (multi_NUBspline_3d_c *spline,
			      double x, double y, double z,
			      complex_float* restrict vals,
			      complex_float* restrict grads);

void
eval_multi_NUBspline_3d_c_vgl (multi_NUBspline_3d_c *spline,
			       double x, double y, double z,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict lapl);

void
eval_multi_NUBspline_3d_c_vgh (multi_NUBspline_3d_c *spline,
			       double x, double y, double z,
			       complex_float* restrict vals,
			       complex_float* restrict grads,
			       complex_float* restrict hess);
#endif
