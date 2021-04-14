/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef NUBSPLINE_EVAL_D_H
#define NUBSPLINE_EVAL_D_H

#include <math.h>
#include <stdio.h>
#include "nubspline_structs.h"

/************************************************************/
/* 1D single-precision, real evaluation functions           */
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
/* 2D single-precision, real evaluation functions           */
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
/* 3D single-precision, real evaluation functions           */
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
