/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_EVAL_D_H
#define MULTI_BSPLINE_EVAL_D_H

#include <math.h>
#include <stdio.h>
#include "multi_bspline_structs.h"

/************************************************************/
/* 1D double-precision, real evaluation functions        */
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
/* 2D double-precision, real evaluation functions        */
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
/* 3D double-precision, real evaluation functions           */
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
