/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_EVAL_Z_H
#define MULTI_BSPLINE_EVAL_Z_H


/************************************************************/
/* 1D double-precision, complex evaluation functions        */
/************************************************************/






/************************************************************/
/* 2D double-precision, complex evaluation functions        */
/************************************************************/




/************************************************************/
/* 3D double-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_UBspline_3d_z (const multi_UBspline_3d_z *spline,
                          double x, double y, double z,
                          complex_double* restrict vals);

void
eval_multi_UBspline_3d_z_vg (const multi_UBspline_3d_z *spline,
                             double x, double y, double z,
                             complex_double* restrict vals,
                             complex_double* restrict grads);

void
eval_multi_UBspline_3d_z_vgl (const multi_UBspline_3d_z *spline,
                              double x, double y, double z,
                              complex_double* restrict vals,
                              complex_double* restrict grads,
                              complex_double* restrict lapl);

void
eval_multi_UBspline_3d_z_vgh (const multi_UBspline_3d_z *spline,
                              double x, double y, double z,
                              complex_double* restrict vals,
                              complex_double* restrict grads,
                              complex_double* restrict hess);



#endif
