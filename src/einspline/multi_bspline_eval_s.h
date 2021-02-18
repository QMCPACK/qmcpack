/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#ifndef MULTI_BSPLINE_EVAL_S_H
#define MULTI_BSPLINE_EVAL_S_H

#include "multi_bspline_structs.h"

/************************************************************/
/* 1D single-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_s (const multi_UBspline_1d_s *spline,
                          float x,
                          float* restrict vals);

void
eval_multi_UBspline_1d_s_vg (const multi_UBspline_1d_s *spline,
                             float x,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_UBspline_1d_s_vgl (const multi_UBspline_1d_s *spline,
                              float x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);


void
eval_multi_UBspline_1d_s_vgh (const multi_UBspline_1d_s *spline,
                              float x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);

/************************************************************/
/* 2D single-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_s(const multi_UBspline_2d_s *spline,
                         double x, double y,
                         float* restrict vals);

void
eval_multi_UBspline_2d_s_vg (const multi_UBspline_2d_s *spline,
                             double x, double y,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_UBspline_2d_s_vgl (const multi_UBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);

void
eval_multi_UBspline_2d_s_vgh (const multi_UBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);

/************************************************************/
/* 3D single-precision, complex evaluation functions        */
/************************************************************/
void
eval_multi_UBspline_3d_s (const multi_UBspline_3d_s *spline,
                          float x, float y, float z,
                          float* restrict vals);

void
eval_multi_UBspline_3d_s_vg (const multi_UBspline_3d_s *spline,
                             float x, float y, float z,
                             float* restrict vals,
                             float* restrict grads);

void
eval_multi_UBspline_3d_s_vgl (const multi_UBspline_3d_s *spline,
                              float x, float y, float z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl);

void
eval_multi_UBspline_3d_s_vgh (const multi_UBspline_3d_s *spline,
                              float x, float y, float z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess);
#endif
