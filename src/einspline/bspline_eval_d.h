/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating b-splines            //
//  copyright (c) 2007 kenneth p. esler, jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_EVAL_D_H
#define BSPLINE_EVAL_D_H

/************************************************************/
/* 1d double-precision, real evaluation functions           */
/************************************************************/

void
eval_UBspline_1d_d (UBspline_1d_d * restrict spline,
                    double x, double* restrict val);

void
eval_UBspline_1d_d_vg (UBspline_1d_d * restrict spline, double x,
                       double* restrict val, double* restrict grad);
void
eval_UBspline_1d_d_vgl (UBspline_1d_d * restrict spline, double x,
                        double* restrict val, double* restrict grad,
                        double* restrict lapl);
void
eval_UBspline_1d_d_vgh (UBspline_1d_d * restrict spline, double x,
                        double* restrict val, double* restrict grad,
                        double* restrict hess);


/************************************************************/
/* 2d double-precision, real evaluation functions           */
/************************************************************/

void
eval_UBspline_2d_d (UBspline_2d_d * restrict spline,
                    double x, double y, double* restrict val);

void
eval_UBspline_2d_d_vg (UBspline_2d_d * restrict spline,
                       double x, double y,
                       double* restrict val, double* restrict grad);

void
eval_UBspline_2d_d_vgl (UBspline_2d_d * restrict spline,
                        double x, double y, double* restrict val,
                        double* restrict grad, double* restrict lapl);

void
eval_UBspline_2d_d_vgh (UBspline_2d_d * restrict spline,
                        double x, double y, double* restrict val,
                        double* restrict grad, double* restrict hess);

/************************************************************/
/* 3d double-precision, real evaluation functions           */
/************************************************************/

void
eval_UBspline_3d_d (UBspline_3d_d * restrict spline,
                    double x, double y, double z,
                    double* restrict val);

void
eval_UBspline_3d_d_vg (UBspline_3d_d * restrict spline,
                       double x, double y, double z,
                       double* restrict val, double* restrict grad);

void
eval_UBspline_3d_d_vgl (UBspline_3d_d * restrict spline,
                        double x, double y, double z,
                        double* restrict val, double* restrict grad, double* restrict lapl);

void
eval_UBspline_3d_d_vgh (UBspline_3d_d * restrict spline,
                        double x, double y, double z,
                        double* restrict val, double* restrict grad, double* restrict hess);

#endif
