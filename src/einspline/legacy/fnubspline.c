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


#include "fnubspline.h"

#include "config.h"
#include "nubspline_create.h"
#include "nubspline.h"

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
  (double *points, int *num_points, NUgrid **grid)
{
  *grid = create_general_grid (points, *num_points);
}

CFUNC void 
F77_FUNC_(fcreate_center_grid,FCREATE_CENTER_GRID)
  (double *start, double *end, double *ratio,
   int *num_points, NUgrid **grid)
{
  *grid = create_center_grid (*start, *end, *ratio, *num_points);
}

CFUNC void
F77_FUNC_(fdestroy_grid,FDESTROY_GRID)
  (NUgrid **grid)
{
  destroy_grid (*grid);
}

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
   float *data, NUBspline_1d_s **spline)
{
  BCtype_s xBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;

  *spline = create_NUBspline_1d_s (*x_grid, xBC, data);
}


CFUNC void
F77_FUNC_(fcreate_nubspline_1d_d,FCREATE_NUBSPLINE_1D_D)
  (NUgrid **x_grid, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   double *data, NUBspline_1d_d **spline)
{
  BCtype_d xBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;

  *spline = create_NUBspline_1d_d (*x_grid, xBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_1d_c,FCREATE_NUBSPLINE_1D_C)
  (NUgrid **x_grid, 
   int *x0_code, complex_float *x0_val, 
   int *x1_code, complex_float *x1_val,
   complex_float *data, NUBspline_1d_c **spline)
{
  BCtype_c xBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = crealf(*x0_val);
  xBC.lVal_i  = cimagf(*x0_val);
  xBC.rVal_r  = crealf(*x1_val);
  xBC.rVal_i  = cimagf(*x1_val);

  *spline = create_NUBspline_1d_c (*x_grid, xBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_1d_z,FCREATE_NUBSPLINE_1D_Z)
  (NUgrid **x_grid, 
   int *x0_code, complex_double *x0_val, 
   int *x1_code, complex_double *x1_val,
   complex_double *data, NUBspline_1d_z **spline)
{
  BCtype_z xBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = creal(*x0_val);
  xBC.lVal_i  = cimag(*x0_val);
  xBC.rVal_r  = creal(*x1_val);
  xBC.rVal_i  = cimag(*x1_val);

  *spline = create_NUBspline_1d_z (*x_grid, xBC, data);
}

////////
// 2D //
////////
CFUNC void
F77_FUNC_(fcreate_nubspline_2d_s,FCREATE_NUBSPLINE_2D_S)
  (NUgrid **x_grid, NUgrid **y_grid, 
   int* x0_code, float *x0_val, int *x1_code, float *x1_val,
   int* y0_code, float *y0_val, int *y1_code, float *y1_val,
   float *data, NUBspline_2d_s **spline)
{
  BCtype_s xBC, yBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;

  *spline = create_NUBspline_2d_s (*x_grid, *y_grid, xBC, yBC, data);
}


CFUNC void
F77_FUNC_(fcreate_nubspline_2d_d,FCREATE_NUBSPLINE_2D_D)
  (NUgrid **x_grid, NUgrid **y_grid,
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   double *data, NUBspline_2d_d **spline)
{
  BCtype_d xBC, yBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;

  *spline = create_NUBspline_2d_d (*x_grid, *y_grid, xBC, yBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_2d_c,FCREATE_NUBSPLINE_2D_C)
  (NUgrid **x_grid, NUgrid **y_grid,
   int *x0_code, complex_float *x0_val, 
   int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, 
   int *y1_code, complex_float *y1_val,
   complex_float *data, NUBspline_2d_c **spline)
{
  BCtype_c xBC, yBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = crealf(*x0_val);
  xBC.lVal_i  = cimagf(*x0_val);
  xBC.rVal_r  = crealf(*x1_val);
  xBC.rVal_i  = cimagf(*x1_val);
  
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal_r  = crealf(*y0_val);
  yBC.lVal_i  = cimagf(*y0_val);
  yBC.rVal_r  = crealf(*y1_val);
  yBC.rVal_i  = cimagf(*y1_val);

  *spline = create_NUBspline_2d_c (*x_grid, *y_grid, xBC, yBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_2d_z,FCREATE_NUBSPLINE_2D_Z)
  (NUgrid **x_grid, NUgrid **y_grid,
   int *x0_code, complex_double *x0_val, 
   int *x1_code, complex_double *x1_val,
   int *y0_code, complex_double *y0_val, 
   int *y1_code, complex_double *y1_val,
   complex_double *data, NUBspline_2d_z **spline)
{
  BCtype_z xBC, yBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = creal(*x0_val);
  xBC.lVal_i  = cimag(*x0_val);
  xBC.rVal_r  = creal(*x1_val);
  xBC.rVal_i  = cimag(*x1_val);

  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal_r  = creal(*y0_val);
  yBC.lVal_i  = cimag(*y0_val);
  yBC.rVal_r  = creal(*y1_val);
  yBC.rVal_i  = cimag(*y1_val);

  *spline = create_NUBspline_2d_z (*x_grid, *y_grid, xBC, yBC, data);
}

////////
// 3D //
////////
CFUNC void
F77_FUNC_(fcreate_nubspline_3d_s,FCREATE_NUBSPLINE_3D_S)
  (NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
   int* x0_code, float *x0_val, int *x1_code, float *x1_val,
   int* y0_code, float *y0_val, int *y1_code, float *y1_val,
   int* z0_code, float *z0_val, int *z1_code, float *z1_val,
   float *data, NUBspline_3d_s **spline)
{
  BCtype_s xBC, yBC, zBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;
  zBC.lCode = (bc_code) *z0_code;
  zBC.rCode = (bc_code) *z1_code;
  zBC.lVal  = *z0_val;
  zBC.rVal  = *z1_val;

  *spline = create_NUBspline_3d_s (*x_grid, *y_grid, *z_grid,
				   xBC, yBC, zBC, data);
}


CFUNC void
F77_FUNC_(fcreate_nubspline_3d_d,FCREATE_NUBSPLINE_3D_D)
  (NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   int* z0_code, float *z0_val, int *z1_code, float *z1_val,
   double *data, NUBspline_3d_d **spline)
{
  BCtype_d xBC, yBC, zBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;
  zBC.lCode = (bc_code) *z0_code;
  zBC.rCode = (bc_code) *z1_code;
  zBC.lVal  = *z0_val;
  zBC.rVal  = *z1_val;

  *spline = create_NUBspline_3d_d (*x_grid, *y_grid, *z_grid,
				   xBC, yBC, zBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_3d_c,FCREATE_NUBSPLINE_3D_C)
  (NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
   int *x0_code, complex_float *x0_val, 
   int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, 
   int *y1_code, complex_float *y1_val,
   int *z0_code, complex_float *z0_val, 
   int *z1_code, complex_float *z1_val,
   complex_float *data, NUBspline_3d_c **spline)
{
  BCtype_c xBC, yBC, zBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = crealf(*x0_val);
  xBC.lVal_i  = cimagf(*x0_val);
  xBC.rVal_r  = crealf(*x1_val);
  xBC.rVal_i  = cimagf(*x1_val);
  
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal_r  = crealf(*y0_val);
  yBC.lVal_i  = cimagf(*y0_val);
  yBC.rVal_r  = crealf(*y1_val);
  yBC.rVal_i  = cimagf(*y1_val);

  zBC.lCode = (bc_code) *z0_code;
  zBC.rCode = (bc_code) *z1_code;
  zBC.lVal_r  = crealf(*z0_val);
  zBC.lVal_i  = cimagf(*z0_val);
  zBC.rVal_r  = crealf(*z1_val);
  zBC.rVal_i  = cimagf(*z1_val);

  *spline = create_NUBspline_3d_c (*x_grid, *y_grid, *z_grid,
				   xBC, yBC, zBC, data);
}

CFUNC void
F77_FUNC_(fcreate_nubspline_3d_z,FCREATE_NUBSPLINE_3D_Z)
  (NUgrid **x_grid, NUgrid **y_grid, NUgrid **z_grid,
   int *x0_code, complex_double *x0_val, 
   int *x1_code, complex_double *x1_val,
   int *y0_code, complex_double *y0_val, 
   int *y1_code, complex_double *y1_val,
   int *z0_code, complex_float *z0_val, 
   int *z1_code, complex_float *z1_val,
   complex_double *data, NUBspline_3d_z **spline)
{
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = creal(*x0_val);
  xBC.lVal_i  = cimag(*x0_val);
  xBC.rVal_r  = creal(*x1_val);
  xBC.rVal_i  = cimag(*x1_val);

  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal_r  = creal(*y0_val);
  yBC.lVal_i  = cimag(*y0_val);
  yBC.rVal_r  = creal(*y1_val);
  yBC.rVal_i  = cimag(*y1_val);

  zBC.lCode = (bc_code) *z0_code;
  zBC.rCode = (bc_code) *z1_code;
  zBC.lVal_r  = creal(*z0_val);
  zBC.lVal_i  = cimag(*z0_val);
  zBC.rVal_r  = creal(*z1_val);
  zBC.rVal_i  = cimag(*z1_val);

  *spline = create_NUBspline_3d_z (*x_grid, *y_grid, *z_grid,
				   xBC, yBC, zBC, data);
}

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
  (NUBspline_1d_s **spline, double *x, float *val)
{
  eval_NUBspline_1d_s (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vg,FEVAL_NUBSPLINE_1D_S_VG)
  (NUBspline_1d_s **spline, double *x, float *val, float *grad)
{
  eval_NUBspline_1d_s_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vgl,FEVAL_NUBSPLINE_1D_S_VGL)
  (NUBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *lapl)
{
  eval_NUBspline_1d_s_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_s_vgh,FEVAL_NUBSPLINE_1D_S_VGH)
  (NUBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *hess)
{
  eval_NUBspline_1d_s_vgh (*spline, *x, val, grad, hess);
}

//////////////////////////////
// 1D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_1d_d,FEVAL_NUBSPLINE_1D_D)
  (NUBspline_1d_d **spline, double *x, double *val)
{
  eval_NUBspline_1d_d (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vg,FEVAL_NUBSPLINE_1D_D_VG)
  (NUBspline_1d_d **spline, double *x, 
   double *val, double *grad)
{
  eval_NUBspline_1d_d_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vgl,FEVAL_NUBSPLINE_1D_D_VGL)
  (NUBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *lapl)
{
  eval_NUBspline_1d_d_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_d_vgh,FEVAL_NUBSPLINE_1D_D_VGH)
  (NUBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *hess)
{
  eval_NUBspline_1d_d_vgh (*spline, *x, val, grad, hess);
}

/////////////////////////////////
// 1D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_1d_c,FEVAL_NUBSPLINE_1D_C)
  (NUBspline_1d_c **spline, double *x, complex_float *val)
{
  eval_NUBspline_1d_c (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vg,FEVAL_NUBSPLINE_1D_C_VG)
  (NUBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad)
{
  eval_NUBspline_1d_c_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vgl,FEVAL_NUBSPLINE_1D_C_VGL)
  (NUBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_NUBspline_1d_c_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_c_vgh,FEVAL_NUBSPLINE_1D_C_VGH)
  (NUBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_NUBspline_1d_c_vgh (*spline, *x, val, grad, hess);
}

/////////////////////////////////
// 1D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nnubspline_1d_z,FEVAL_NNUBSPLINE_1D_Z)
  (NUBspline_1d_z **spline, double *x, complex_double *val)
{
  eval_NUBspline_1d_z (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vg,FEVAL_NUBSPLINE_1D_Z_VG)
  (NUBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad)
{
  eval_NUBspline_1d_z_vg (*spline, *x, val, grad);
}
 
CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vgl,FEVAL_NUBSPLINE_1D_Z_VGL)
  (NUBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_NUBspline_1d_z_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_1d_z_vgh,FEVAL_NUBSPLINE_1D_Z_VGH)
  (NUBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_NUBspline_1d_z_vgh (*spline, *x, val, grad, hess);
}

//////////////////////////////
// 2D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_s,FEVAL_NUBSPLINE_2D_S)
  (NUBspline_2d_s **spline, double *x, double *y, float *val)
{
  eval_NUBspline_2d_s (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vg,FEVAL_NUBSPLINE_2D_S_VG)
  (NUBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad)
{
  eval_NUBspline_2d_s_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vgl,FEVAL_NUBSPLINE_2D_S_VGL)
  (NUBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float* lapl)
{
  eval_NUBspline_2d_s_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_s_vgh,FEVAL_NUBSPLINE_2D_S_VGH)
  (NUBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float *hess)
{
  eval_NUBspline_2d_s_vgh (*spline, *x, *y, val, grad, hess);
}

//////////////////////////////
// 2D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_d,FEVAL_NUBSPLINE_2D_D)
  (NUBspline_2d_d **spline, double *x, double *y, double *val)
{
  eval_NUBspline_2d_d (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vg,FEVAL_NUBSPLINE_2D_D_VG)
  (NUBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad)
{
  eval_NUBspline_2d_d_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vgl,FEVAL_NUBSPLINE_2D_D_VGL)
  (NUBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *lapl)
{
  eval_NUBspline_2d_d_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_d_vgh,FEVAL_NUBSPLINE_2D_D_VGH)
  (NUBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *hess)
{
  eval_NUBspline_2d_d_vgl (*spline, *x, *y, val, grad, hess);
}

/////////////////////////////////
// 2D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_c,FEVAL_NUBSPLINE_2D_C)
  (NUBspline_2d_c **spline, double *x, double *y, complex_float *val)
{
  eval_NUBspline_2d_c (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vg,FEVAL_NUBSPLINE_2D_C_VG)
  (NUBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad)
{
  eval_NUBspline_2d_c_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vgl,FEVAL_NUBSPLINE_2D_C_VGL)
  (NUBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_NUBspline_2d_c_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_c_vgh,FEVAL_NUBSPLINE_2D_C_VGH)
  (NUBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_NUBspline_2d_c_vgh (*spline, *x, *y, val, grad, hess);
}

/////////////////////////////////
// 2D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_2d_z,FEVAL_NUBSPLINE_2D_Z)
  (NUBspline_2d_z **spline, double *x, double *y, complex_double *val)
{
  eval_NUBspline_2d_z (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vg,FEVAL_NUBSPLINE_2D_Z_VG)
  (NUBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad)
{
  eval_NUBspline_2d_z_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vgl,FEVAL_NUBSPLINE_2D_Z_VGL)
  (NUBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_NUBspline_2d_z_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_2d_z_vgh,FEVAL_NUBSPLINE_2D_Z_VGH)
  (NUBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_NUBspline_2d_z_vgh (*spline, *x, *y, val, grad, hess);
}



//////////////////////////////
// 3D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_s,FEVAL_NUBSPLINE_3D_S)
  (NUBspline_3d_s **spline, double *x, double *y, double *z,
   float *val)
{
  eval_NUBspline_3d_s (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vg,FEVAL_NUBSPLINE_3D_S_VG)
  (NUBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad)
{
  eval_NUBspline_3d_s_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vgl,FEVAL_NUBSPLINE_3D_S_VGL)
  (NUBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad, float* lapl)
{
  eval_NUBspline_3d_s_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_s_vgh,FEVAL_NUBSPLINE_3D_S_VGH)
  (NUBspline_3d_s **spline, double *x, double *y, double *z, 
   float *val, float *grad, float *hess)
{
  eval_NUBspline_3d_s_vgh (*spline, *x, *y, *z, val, grad, hess);
}

//////////////////////////////
// 3D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_d,FEVAL_NUBSPLINE_3D_D)
  (NUBspline_3d_d **spline, double *x, double *y, double *z,
   double *val)
{
  eval_NUBspline_3d_d (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vg,FEVAL_NUBSPLINE_3D_D_VG)
  (NUBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad)
{
  eval_NUBspline_3d_d_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vgl,FEVAL_NUBSPLINE_3D_D_VGL)
  (NUBspline_3d_d **spline, double *x, double *y, double *z,  
   double *val, double *grad, double *lapl)
{
  eval_NUBspline_3d_d_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_d_vgh,FEVAL_NUBSPLINE_3D_D_VGH)
  (NUBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad, double *hess)
{
  eval_NUBspline_3d_d_vgl (*spline, *x, *y, *z, val, grad, hess);
}

/////////////////////////////////
// 3D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_c,FEVAL_NUBSPLINE_3D_C)
  (NUBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val)
{
  eval_NUBspline_3d_c (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vg,FEVAL_NUBSPLINE_3D_C_VG)
  (NUBspline_3d_c **spline, double *x, double *y, double *z, 
   complex_float *val, complex_float *grad)
{
  eval_NUBspline_3d_c_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vgl,FEVAL_NUBSPLINE_3D_C_VGL)
  (NUBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_NUBspline_3d_c_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_c_vgh,FEVAL_NUBSPLINE_3D_C_VGH)
  (NUBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_NUBspline_3d_c_vgh (*spline, *x, *y, *z, val, grad, hess);
}

/////////////////////////////////
// 3D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_nubspline_3d_z,FEVAL_NUBSPLINE_3D_Z)
  (NUBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val)
{
  eval_NUBspline_3d_z (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vg,FEVAL_NUBSPLINE_3D_Z_VG)
  (NUBspline_3d_z **spline, double *x, double *y, double *z, 
   complex_double *val, complex_double *grad)
{
  eval_NUBspline_3d_z_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vgl,FEVAL_NUBSPLINE_3D_Z_VGL)
  (NUBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_NUBspline_3d_z_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_nubspline_3d_z_vgh,FEVAL_NUBSPLINE_3D_Z_VGH)
  (NUBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_NUBspline_3d_z_vgh (*spline, *x, *y, *z, val, grad, hess);
}

