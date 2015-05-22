/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating b-splines            //
//  copyright (c) 2007 kenneth p. esler, jr.                               //
//                                                                         //
//  this program is free software; you can redistribute it and/or modify   //
//  it under the terms of the gnu general public license as published by   //
//  the free software foundation; either version 2 of the license, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  this program is distributed in the hope that it will be useful,        //
//  but without any warranty; without even the implied warranty of         //
//  merchantability or fitness for a particular purpose.  see the          //
//  gnu general public license for more details.                           //
//                                                                         //
//  you should have received a copy of the gnu general public license      //
//  along with this program; if not, write to the free software            //
//  foundation, inc., 51 franklin street, fifth floor,                     //
//  boston, ma  02110-1301  usa                                            //
/////////////////////////////////////////////////////////////////////////////

#ifndef BSPLINE_EVAL_D_H
#define BSPLINE_EVAL_D_H

/************************************************************/
/* 1d double-precision, real evaulation functions           */
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
/* 2d double-precision, real evaulation functions           */
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
/* 3d double-precision, real evaulation functions           */
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
