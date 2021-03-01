/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_NUBSPLINE_EVAL_Z_H
#define MULTI_NUBSPLINE_EVAL_Z_H


/************************************************************/
/* 1D double-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_NUBspline_1d_z (multi_NUBspline_1d_z *spline,
                           double x,
                           complex_double* restrict vals);

void
eval_multi_NUBspline_1d_z_vg (multi_NUBspline_1d_z *spline,
                              double x,
                              complex_double* restrict vals,
                              complex_double* restrict grads);

void
eval_multi_NUBspline_1d_z_vgl (multi_NUBspline_1d_z *spline,
                               double x,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict lapl);


void
eval_multi_NUBspline_1d_z_vgh (multi_NUBspline_1d_z *spline,
                               double x,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict hess);


/************************************************************/
/* 2D double-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_NUBspline_2d_z (multi_NUBspline_2d_z *spline,
                           double x, double y,
                           complex_double* restrict vals);

void
eval_multi_NUBspline_2d_z_vg (multi_NUBspline_2d_z *spline,
                              double x, double y,
                              complex_double* restrict vals,
                              complex_double* restrict grads);

void
eval_multi_NUBspline_2d_z_vgl (multi_NUBspline_2d_z *spline,
                               double x, double y,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict lapl);

void
eval_multi_NUBspline_2d_z_vgh (multi_NUBspline_2d_z *spline,
                               double x, double y,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict hess);

/************************************************************/
/* 3D double-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_NUBspline_3d_z (multi_NUBspline_3d_z *spline,
                           double x, double y, double z,
                           complex_double* restrict vals);

void
eval_multi_NUBspline_3d_z_vg (multi_NUBspline_3d_z *spline,
                              double x, double y, double z,
                              complex_double* restrict vals,
                              complex_double* restrict grads);

void
eval_multi_NUBspline_3d_z_vgl (multi_NUBspline_3d_z *spline,
                               double x, double y, double z,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict lapl);

void
eval_multi_NUBspline_3d_z_vgh (multi_NUBspline_3d_z *spline,
                               double x, double y, double z,
                               complex_double* restrict vals,
                               complex_double* restrict grads,
                               complex_double* restrict hess);


#endif
