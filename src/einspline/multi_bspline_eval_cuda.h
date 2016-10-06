//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
// 		      Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef MULTI_BSPLINE_EVAL_CUDA_H
#define MULTI_BSPLINE_EVAL_CUDA_H

#include "multi_bspline_structs_cuda.h"


////////
// 1D //
////////

// Single-precision real
extern "C" void
eval_multi_multi_UBspline_1d_s_cuda
(const multi_UBspline_1d_s_cuda *spline, float *pos_d, float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_1d_s_vgl_cuda
(const multi_UBspline_1d_s_cuda *spline, float *pos_d,
 float *vals_d[], float *grads_d[], float *lapl_d[], int num);

// Double-precision real
extern "C" void
eval_multi_multi_UBspline_1d_d_cuda
(const multi_UBspline_1d_d_cuda *spline, double *pos_d, double *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_1d_d_vgl_cuda
(const multi_UBspline_1d_d_cuda *spline, double *pos_d,
 double *vals_d[], double *grad_lapl_d[], int num, int row_stride);


// Single-precision complex
extern "C" void
eval_multi_multi_UBspline_1d_c_cuda
(const multi_UBspline_1d_c_cuda *spline,
 float *pos_d, complex_float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_1d_c_vgl_cuda
(const multi_UBspline_1d_c_cuda *spline, float *pos_d,
 complex_float *vals_d[], complex_float *grad_lapl_d[], int num, int row_stride);


// Doublele-precision complex
extern "C" void
eval_multi_multi_UBspline_1d_z_cuda
(const multi_UBspline_1d_z_cuda *spline,
 double *pos_d, complex_double *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_1d_z_vgl_cuda
(const multi_UBspline_1d_z_cuda *spline, double *pos_d,
 complex_double *vals_d[], complex_double *grad_lapl_d[], int num, int row_stride);




////////
// 3D //
////////
// Single-precision real
extern "C" void
eval_multi_multi_UBspline_3d_s_cuda
(const multi_UBspline_3d_s_cuda *spline, float *pos_d, float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_s_sign_cuda
(const multi_UBspline_3d_s_cuda *spline, float *pos_d, float *sign_d,
 float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_s_vgh_cuda
(const multi_UBspline_3d_s_cuda *spline,
 float *pos_d, float *vals_d[], float *grads_d[], float *hess_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_s_vgl_cuda
(const multi_UBspline_3d_s_cuda *spline, float *pos_d, float *Linv_d,
 float *vals_d[], float *grad_lapl_d[], int num, int row_stride);

extern "C" void
eval_multi_multi_UBspline_3d_s_vgl_sign_cuda
(const multi_UBspline_3d_s_cuda *spline, float *pos_d, float *sign_d, float *Linv_d,
 float *vals_d[], float *grad_lapl_d[], int num, int row_stride);



// Double-precision real
extern "C" void
eval_multi_multi_UBspline_3d_d_cuda
(const multi_UBspline_3d_d_cuda *spline, double *pos_d, double *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_d_sign_cuda
(const multi_UBspline_3d_d_cuda *spline, double *pos_d, double *sign_d,
 double *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_d_vgh_cuda
(const multi_UBspline_3d_d_cuda *spline,
 double *pos_d, double *vals_d[], double *grads_d[], double *hess_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_d_vgl_cuda
(const multi_UBspline_3d_d_cuda *spline, double *pos_d, double *Linv_d,
 double *vals_d[], double *grad_lapl_d[], int num, int row_stride);

extern "C" void
eval_multi_multi_UBspline_3d_d_vgl_sign_cuda
(const multi_UBspline_3d_d_cuda *spline, double *pos_d, double *sign_d, double *Linv_d,
 double *vals_d[], double *grad_lapl_d[], int num, int row_stride);


// Single-precision complex
extern "C" void
eval_multi_multi_UBspline_3d_c_cuda
(const multi_UBspline_3d_c_cuda *spline,
 float *pos_d, complex_float *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_c_vgh_cuda
(const multi_UBspline_3d_c_cuda *spline, float *pos_d,
 complex_float *vals_d[], complex_float *grads_d[],
 complex_float *hess_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_c_vgl_cuda
(const multi_UBspline_3d_c_cuda *spline, float *pos_d, float *Linv_d,
 complex_float *vals_d[], complex_float *grad_lapl_d[], int num, int row_stride);


// Doublele-precision complex
extern "C" void
eval_multi_multi_UBspline_3d_z_cuda
(const multi_UBspline_3d_z_cuda *spline,
 double *pos_d, complex_double *vals_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_z_vgh_cuda
(const multi_UBspline_3d_z_cuda *spline, double *pos_d,
 complex_double *vals_d[], complex_double *grads_d[],
 complex_double *hess_d[], int num);

extern "C" void
eval_multi_multi_UBspline_3d_z_vgl_cuda
(const multi_UBspline_3d_z_cuda *spline, double *pos_d, double *Linv_d,
 complex_double *vals_d[], complex_double *grad_lapl_d[], int num, int row_stride);


#endif
