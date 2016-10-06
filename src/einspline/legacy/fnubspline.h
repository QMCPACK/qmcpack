//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign 
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef F_NUBSPLINE_H
#define F_NUBSPLINE_H

#include "config.h"
#include "nugrid.h"
#include "nubspline_structs.h"

#ifdef __cplusplus
#define CFUNC extern "C" /* Avoid name mangling in C++ */
#else
#define CFUNC
#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////                    Grid Creation routines                    ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

CFUNC void
F77_FUNC_(fcreate_general_grid,FCREATE_GENERAL_GRID)
(double *points, int *num_points, NUgrid **grid);

CFUNC void
F77_FUNC_(fcreate_center_grid,FCREATE_CENTER_GRID)
(double *start, double *end, double *ratio,
 int *num_points, NUgrid **grid);

CFUNC void
F77_FUNC_(fdestroy_grid,FDESTROY_GRID)
(NUgrid **grid);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////            Nonuniform spline creation routines               ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

////////
// 1D //
////////
CFUNC void
F77_FUNC_(fcreate_nubspline_1d_s,FCREATE_NUBSPLINE_1D_S)
(NUgrid **x_grid,
 int* x0_code, float *x0_val, int *x1_code, float *x1_val,
 float *data, NUBspline_1d_s **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_1d_d,FCREATE_NUBSPLINE_1D_D)
(NUgrid **x_grid,
 int *x0_code, double *x0_val, int *x1_code, double *x1_val,
 double *data, NUBspline_1d_d **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_1d_c,FCREATE_NUBSPLINE_1D_C)
(NUgrid **x_grid,
 int *x0_code, complex_float *x0_val,
 int *x1_code, complex_float *x1_val,
 complex_float *data, NUBspline_1d_c **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_1d_z,FCREATE_NUBSPLINE_1D_Z)
(NUgrid **x_grid,
 int *x0_code, complex_double *x0_val,
 int *x1_code, complex_double *x1_val,
 complex_double *data, NUBspline_1d_z **spline);

////////
// 2D //
////////
CFUNC void
F77_FUNC_(fcreate_nubspline_2d_s,FCREATE_NUBSPLINE_2D_S)
(NUgrid **x_grid, NUgrid **y_grid,
 int* x0_code, float *x0_val, int *x1_code, float *x1_val,
 int* y0_code, float *y0_val, int *y1_code, float *y1_val,
 float *data, NUBspline_2d_s **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_2d_d,FCREATE_NUBSPLINE_2D_D)
(NUgrid **x_grid, NUgrid **y_grid,
 int *x0_code, double *x0_val, int *x1_code, double *x1_val,
 int *y0_code, double *y0_val, int *y1_code, double *y1_val,
 double *data, NUBspline_2d_d **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_2d_c,FCREATE_NUBSPLINE_2D_C)
(NUgrid **x_grid, NUgrid **y_grid,
 int *x0_code, complex_float *x0_val,
 int *x1_code, complex_float *x1_val,
 int *y0_code, complex_float *y0_val,
 int *y1_code, complex_float *y1_val,
 complex_float *data, NUBspline_2d_c **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_2d_z,FCREATE_NUBSPLINE_2D_Z)
(NUgrid **x_grid, NUgrid **y_grid,
 int *x0_code, complex_double *x0_val,
 int *x1_code, complex_double *x1_val,
 int *y0_code, complex_double *y0_val,
 int *y1_code, complex_double *y1_val,
 complex_double *data, NUBspline_2d_z **spline);

////////
// 3D //
////////
CFUNC void
F77_FUNC_(fcreate_nubspline_3d_s,FCREATE_NUBSPLINE_3D_S)
(NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
 int* x0_code, float *x0_val, int *x1_code, float *x1_val,
 int* y0_code, float *y0_val, int *y1_code, float *y1_val,
 int* z0_code, float *z0_val, int *z1_code, float *z1_val,
 float *data, NUBspline_3d_s **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_3d_d,FCREATE_NUBSPLINE_3D_D)
(NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
 int *x0_code, double *x0_val, int *x1_code, double *x1_val,
 int *y0_code, double *y0_val, int *y1_code, double *y1_val,
 int* z0_code, float *z0_val, int *z1_code, float *z1_val,
 double *data, NUBspline_3d_d **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_3d_c,FCREATE_NUBSPLINE_3D_C)
(NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
 int *x0_code, complex_float *x0_val,
 int *x1_code, complex_float *x1_val,
 int *y0_code, complex_float *y0_val,
 int *y1_code, complex_float *y1_val,
 int *z0_code, complex_float *z0_val,
 int *z1_code, complex_float *z1_val,
 complex_float *data, NUBspline_3d_c **spline);

CFUNC void
F77_FUNC_(fcreate_nubspline_3d_z,FCREATE_NUBSPLINE_3D_Z)
(NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
 int *x0_code, complex_double *x0_val,
 int *x1_code, complex_double *x1_val,
 int *y0_code, complex_double *y0_val,
 int *y1_code, complex_double *y1_val,
 int *z0_code, complex_float *z0_val,
 int *z1_code, complex_float *z1_val,
 complex_double *data, NUBspline_3d_z **spline);


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////           Nonuniform spline evaluation routines              ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////
// 1D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_1d_s,FEVAL_NUBSPLINE_1D_S)
(NUBspline_1d_s **spline, double *x, float *val);

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vg,FEVAL_NUBSPLINE_1D_S_VG)
(NUBspline_1d_s **spline, double *x, float *val, float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vgl,FEVAL_NUBSPLINE_1D_S_VGL)
(NUBspline_1d_s **spline, double *x,
 float *val, float *grad, float *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vgh,FEVAL_NUBSPLINE_1D_S_VGH)
(NUBspline_1d_s **spline, double *x,
 float *val, float *grad, float *hess);

//////////////////////////////
// 1D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_1d_d,FEVAL_NUBSPLINE_1D_D)
(NUBspline_1d_d **spline, double *x, double *val);

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vg,FEVAL_NUBSPLINE_1D_D_VG)
(NUBspline_1d_d **spline, double *x,
 double *val, double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vgl,FEVAL_NUBSPLINE_1D_D_VGL)
(NUBspline_1d_d **spline, double *x,
 double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vgh,FEVAL_NUBSPLINE_1D_D_VGH)
(NUBspline_1d_d **spline, double *x,
 double *val, double *grad, double *hess);

/////////////////////////////////
// 1D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_1d_c,FEVAL_NUBSPLINE_1D_C)
(NUBspline_1d_c **spline, double *x, complex_float *val);

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vg,FEVAL_NUBSPLINE_1D_C_VG)
(NUBspline_1d_c **spline, double *x,
 complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vgl,FEVAL_NUBSPLINE_1D_C_VGL)
(NUBspline_1d_c **spline, double *x,
 complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vgh,FEVAL_NUBSPLINE_1D_C_VGH)
(NUBspline_1d_c **spline, double *x,
 complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 1D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nnubspline_1d_z,FEVAL_NNUBSPLINE_1D_Z)
(NUBspline_1d_z **spline, double *x, complex_double *val);

CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vg,FEVAL_NUBSPLINE_1D_Z_VG)
(NUBspline_1d_z **spline, double *x,
 complex_double *val, complex_double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vgl,FEVAL_NUBSPLINE_1D_Z_VGL)
(NUBspline_1d_z **spline, double *x,
 complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vgh,FEVAL_NUBSPLINE_1D_Z_VGH)
(NUBspline_1d_z **spline, double *x,
 complex_double *val, complex_double *grad, complex_double *hess);

//////////////////////////////
// 2D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_s,FEVAL_NUBSPLINE_2D_S)
(NUBspline_2d_s **spline, double *x, double *y, float *val);

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vg,FEVAL_NUBSPLINE_2D_S_VG)
(NUBspline_2d_s **spline, double *x, double *y,
 float *val, float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vgl,FEVAL_NUBSPLINE_2D_S_VGL)
(NUBspline_2d_s **spline, double *x, double *y,
 float *val, float *grad, float* lapl);

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vgh,FEVAL_NUBSPLINE_2D_S_VGH)
(NUBspline_2d_s **spline, double *x, double *y,
 float *val, float *grad, float *hess);

//////////////////////////////
// 2D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_d,FEVAL_NUBSPLINE_2D_D)
(NUBspline_2d_d **spline, double *x, double *y, double *val);

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vg,FEVAL_NUBSPLINE_2D_D_VG)
(NUBspline_2d_d **spline, double *x, double *y,
 double *val, double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vgl,FEVAL_NUBSPLINE_2D_D_VGL)
(NUBspline_2d_d **spline, double *x, double *y,
 double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vgh,FEVAL_NUBSPLINE_2D_D_VGH)
(NUBspline_2d_d **spline, double *x, double *y,
 double *val, double *grad, double *hess);

/////////////////////////////////
// 2D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_c,FEVAL_NUBSPLINE_2D_C)
(NUBspline_2d_c **spline, double *x, double *y, complex_float *val);

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vg,FEVAL_NUBSPLINE_2D_C_VG)
(NUBspline_2d_c **spline, double *x, double *y,
 complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vgl,FEVAL_NUBSPLINE_2D_C_VGL)
(NUBspline_2d_c **spline, double *x, double *y,
 complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vgh,FEVAL_NUBSPLINE_2D_C_VGH)
(NUBspline_2d_c **spline, double *x, double *y,
 complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 2D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_z,FEVAL_NUBSPLINE_2D_Z)
(NUBspline_2d_z **spline, double *x, double *y, complex_double *val);

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vg,FEVAL_NUBSPLINE_2D_Z_VG)
(NUBspline_2d_z **spline, double *x, double *y,
 complex_double *val, complex_double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vgl,FEVAL_NUBSPLINE_2D_Z_VGL)
(NUBspline_2d_z **spline, double *x, double *y,
 complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vgh,FEVAL_NUBSPLINE_2D_Z_VGH)
(NUBspline_2d_z **spline, double *x, double *y,
 complex_double *val, complex_double *grad, complex_double *hess);

//////////////////////////////
// 3D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_s,FEVAL_NUBSPLINE_3D_S)
(NUBspline_3d_s **spline, double *x, double *y, double *z,
 float *val);

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vg,FEVAL_NUBSPLINE_3D_S_VG)
(NUBspline_3d_s **spline, double *x, double *y, double *z,
 float *val, float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vgl,FEVAL_NUBSPLINE_3D_S_VGL)
(NUBspline_3d_s **spline, double *x, double *y, double *z,
 float *val, float *grad, float* lapl);

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vgh,FEVAL_NUBSPLINE_3D_S_VGH)
(NUBspline_3d_s **spline, double *x, double *y, double *z,
 float *val, float *grad, float *hess);

//////////////////////////////
// 3D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_d,FEVAL_NUBSPLINE_3D_D)
(NUBspline_3d_d **spline, double *x, double *y, double *z,
 double *val);

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vg,FEVAL_NUBSPLINE_3D_D_VG)
(NUBspline_3d_d **spline, double *x, double *y, double *z,
 double *val, double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vgl,FEVAL_NUBSPLINE_3D_D_VGL)
(NUBspline_3d_d **spline, double *x, double *y, double *z,
 double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vgh,FEVAL_NUBSPLINE_3D_D_VGH)
(NUBspline_3d_d **spline, double *x, double *y, double *z,
 double *val, double *grad, double *hess);

/////////////////////////////////
// 3D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_c,FEVAL_NUBSPLINE_3D_C)
(NUBspline_3d_c **spline, double *x, double *y, double *z,
 complex_float *val);

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vg,FEVAL_NUBSPLINE_3D_C_VG)
(NUBspline_3d_c **spline, double *x, double *y, double *z,
 complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vgl,FEVAL_NUBSPLINE_3D_C_VGL)
(NUBspline_3d_c **spline, double *x, double *y, double *z,
 complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vgh,FEVAL_NUBSPLINE_3D_C_VGH)
(NUBspline_3d_c **spline, double *x, double *y, double *z,
 complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 3D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_z,FEVAL_NUBSPLINE_3D_Z)
(NUBspline_3d_z **spline, double *x, double *y, double *z,
 complex_double *val);

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vg,FEVAL_NUBSPLINE_3D_Z_VG)
(NUBspline_3d_z **spline, double *x, double *y, double *z,
 complex_double *val, complex_double *grad);

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vgl,FEVAL_NUBSPLINE_3D_Z_VGL)
(NUBspline_3d_z **spline, double *x, double *y, double *z,
 complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vgh,FEVAL_NUBSPLINE_3D_Z_VGH)
(NUBspline_3d_z **spline, double *x, double *y, double *z,
 complex_double *val, complex_double *grad, complex_double *hess);

#endif
