#ifndef FMULTI_BSPLINE_H
#define FMULTI_BSPLINE_H

#include "config.h"
#include "bspline_base.h"
#include "bspline_create.h"

#ifdef __cplusplus
#define CFUNC extern "C" /* Avoid name mangling in C++ */
#else
#define CFUNC
#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////                       Creation routines                      ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

////////
// 1D //
////////
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_s,FCREATE_MULTI_UBSPLINE_1D_S)
  (double   *x0, double    *x1, int   *num_x, 
   int *x0_code, float *x0_val, int *x1_code, float *x1_val,
   int *num_spline, multi_UBspline_1d_s **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_d,FCREATE_MULTI_UBSPLINE_1D_D)
  (double   *x0, double     *x1, int   *num_x, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *num_splines, multi_UBspline_1d_d **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_c,FCREATE_MULTI_UBSPLINE_1D_C)
  (double   *x0, double    *x1, int   *num_x, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *num_splines, multi_UBspline_1d_c **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_z,FCREATE_MULTI_UBSPLINE_1D_Z)
  (double   *x0, double     *x1, int   *num_x, 
   int *x0_code, complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *num_splines, multi_UBspline_1d_z **spline);

CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_s,FSET_MULTI_UBSPLINE_1D_S)
  (multi_UBspline_1d_s **spline, int *spline_num, float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_d,FSET_MULTI_UBSPLINE_1D_D)
  (multi_UBspline_1d_d **spline, int *spline_num, double *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_c,FSET_MULTI_UBSPLINE_1D_C)
  (multi_UBspline_1d_c **spline, int *spline_num, complex_float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_z,FSET_MULTI_UBSPLINE_1D_Z)
  (multi_UBspline_1d_z **spline, int *spline_num, complex_double *data);

////////
// 2D //
////////
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_s,FCREATE_MULTI_UBSPLINE_2D_S)
  (double   *x0, double    *x1, int   *num_x, 
   double   *y0, double    *y1, int   *num_y, 
   int *x0_code, float *x0_val, int *x1_code, float *x1_val,
   int *y0_code, float *y0_val, int *y1_code, float *y1_val,
   int *num_splines, multi_UBspline_2d_s **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_d,FCREATE_MULTI_UBSPLINE_2D_D)
  (double   *x0, double     *x1, int   *num_x, 
   double   *y0, double     *y1, int   *num_y, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   int *num_splines, multi_UBspline_2d_d **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_c,FCREATE_MULTI_UBSPLINE_2D_C)
  (double   *x0, double    *x1, int   *num_x, 
   double   *y0, double    *y1, int   *num_y, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, int *y1_code, complex_float *y1_val,
   int *num_splines, multi_UBspline_2d_c **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_z,FCREATE_MULTI_UBSPLINE_2D_Z)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   int *x0_code, complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *y0_code, complex_double *y0_val, int *y1_code, complex_double *y1_val,
   int *num_splines, multi_UBspline_2d_z **spline);

CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_s,FSET_MULTI_UBSPLINE_2D_S)
  (multi_UBspline_2d_s **spline, int *spline_num, float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_d,FSET_MULTI_UBSPLINE_2D_D)
  (multi_UBspline_2d_d **spline, int *spline_num, double *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_c,FSET_MULTI_UBSPLINE_2D_C)
  (multi_UBspline_2d_c **spline, int *spline_num, complex_float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_z,FSET_MULTI_UBSPLINE_2D_Z)
  (multi_UBspline_2d_z **spline, int *spline_num, complex_double *data);

////////
// 3D //
////////
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_s,FCREATE_MULTI_UBSPLINE_3D_S)
  (double   *x0, double    *x1, int   *num_x, 
   double   *y0, double    *y1, int   *num_y, 
   double   *z0, double    *z1, int   *num_z, 
   int *x0_code, float *x0_val, int *x1_code, float *x1_val,
   int *y0_code, float *y0_val, int *y1_code, float *y1_val,
   int *z0_code, float *z0_val, int *z1_code, float *z1_val,
   int *num_splines, multi_UBspline_3d_s **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_d,FCREATE_MULTI_UBSPLINE_3D_D)
  (double   *x0, double     *x1, int   *num_x, 
   double   *y0, double     *y1, int   *num_y, 
   double   *z0, double     *z1, int   *num_z, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   int *z0_code, double *z0_val, int *z1_code, double *z1_val,
   int *num_splines, multi_UBspline_3d_d **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_c,FCREATE_MULTI_UBSPLINE_3D_C)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   double *z0, double *z1, int *num_z, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, int *y1_code, complex_float *y1_val,
   int *z0_code, complex_float *z0_val, int *z1_code, complex_float *z1_val,
   int *num_splines, multi_UBspline_3d_c **spline);
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_z,FCREATE_MULTI_UBSPLINE_3D_Z)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   double *z0, double *z1, int *num_z, 
   int *x0_code,  complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *y0_code,  complex_double *y0_val, int *y1_code, complex_double *y1_val,
   int *z0_code,  complex_double *z0_val, int *z1_code, complex_double *z1_val,
   int *num_splines, multi_UBspline_3d_z **spline);

CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_s,FSET_MULTI_UBSPLINE_3D_S)
  (multi_UBspline_3d_s **spline, int *spline_num, float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_d,FSET_MULTI_UBSPLINE_3D_D)
  (multi_UBspline_3d_d **spline, int *spline_num, double *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_c,FSET_MULTI_UBSPLINE_3D_C)
  (multi_UBspline_3d_c **spline, int *spline_num, complex_float *data);
CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_z,FSET_MULTI_UBSPLINE_3D_Z)
  (multi_UBspline_3d_z **spline, int *spline_num, complex_double *data);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////                      Destruction routine                     ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
CFUNC void
F77_FUNC_(fdestroy_bspline,FDESTROY_BSPLINE)
  (Bspline **spline);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
////                      Evaluation routines                     ////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////
// 1D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_ubspline_1d_s,FEVAL_UBSPLINE_1D_S)
  (multi_UBspline_1d_s **spline, double *x, float *val);

CFUNC void
F77_FUNC_(feval_ubspline_1d_s_vg,FEVAL_UBSPLINE_1D_S_VG)
  (multi_UBspline_1d_s **spline, double *x, float *val, float *grad);

CFUNC void
F77_FUNC_(feval_ubspline_1d_s_vgl,FEVAL_UBSPLINE_1D_S_VGL)
  (multi_UBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *lapl);

CFUNC void
F77_FUNC_(feval_ubspline_1d_s_vgh,FEVAL_UBSPLINE_1D_S_VGH)
  (multi_UBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *hess);

//////////////////////////////
// 1D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d,FEVAL_MULTI_UBSPLINE_1D_D)
  (multi_UBspline_1d_d **spline, double *x, double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vg,FEVAL_MULTI_UBSPLINE_1D_D_VG)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vgl,FEVAL_MULTI_UBSPLINE_1D_D_VGL)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vgh,FEVAL_MULTI_UBSPLINE_1D_D_VGH)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *hess);

/////////////////////////////////
// 1D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c,FEVAL_MULTI_UBSPLINE_1D_C)
  (multi_UBspline_1d_c **spline, double *x, complex_float *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vg,FEVAL_MULTI_UBSPLINE_1D_C_VG)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vgl,FEVAL_MULTI_UBSPLINE_1D_C_VGL)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vgh,FEVAL_MULTI_UBSPLINE_1D_C_VGH)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 1D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z,FEVAL_MULTI_UBSPLINE_1D_Z)
  (multi_UBspline_1d_z **spline, double *x, complex_double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vg,FEVAL_MULTI_UBSPLINE_1D_Z_VG)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad);
 
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vgl,FEVAL_MULTI_UBSPLINE_1D_Z_VGL)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vgh,FEVAL_MULTI_UBSPLINE_1D_Z_VGH)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *hess);

//////////////////////////////
// 2D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s,FEVAL_MULTI_UBSPLINE_2D_S)
  (multi_UBspline_2d_s **spline, double *x, double *y, float *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vg,FEVAL_MULTI_UBSPLINE_2D_S_VG)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vgl,FEVAL_MULTI_UBSPLINE_2D_S_VGL)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float* lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vgh,FEVAL_MULTI_UBSPLINE_2D_S_VGH)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float *hess);

//////////////////////////////
// 2D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d,FEVAL_MULTI_UBSPLINE_2D_D)
  (multi_UBspline_2d_d **spline, double *x, double *y, double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vg,FEVAL_MULTI_UBSPLINE_2D_D_VG)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vgl,FEVAL_MULTI_UBSPLINE_2D_D_VGL)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vgh,FEVAL_MULTI_UBSPLINE_2D_D_VGH)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *hess);

/////////////////////////////////
// 2D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c,FEVAL_MULTI_UBSPLINE_2D_C)
  (multi_UBspline_2d_c **spline, double *x, double *y, complex_float *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vg,FEVAL_MULTI_UBSPLINE_2D_C_VG)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vgl,FEVAL_MULTI_UBSPLINE_2D_C_VGL)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vgh,FEVAL_MULTI_UBSPLINE_2D_C_VGH)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 2D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z,FEVAL_MULTI_UBSPLINE_2D_Z)
  (multi_UBspline_2d_z **spline, double *x, double *y, complex_double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vg,FEVAL_MULTI_UBSPLINE_2D_Z_VG)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vgl,FEVAL_MULTI_UBSPLINE_2D_Z_VGL)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vgh,FEVAL_MULTI_UBSPLINE_2D_Z_VGH)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *hess);


//////////////////////////////
// 3D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s,FEVAL_MULTI_UBSPLINE_3D_S)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vg,FEVAL_MULTI_UBSPLINE_3D_S_VG)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vgl,FEVAL_MULTI_UBSPLINE_3D_S_VGL)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad, float* lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vgh,FEVAL_MULTI_UBSPLINE_3D_S_VGH)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z, 
   float *val, float *grad, float *hess);

//////////////////////////////
// 3D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d,FEVAL_MULTI_UBSPLINE_3D_D)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vg,FEVAL_MULTI_UBSPLINE_3D_D_VG)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vgl,FEVAL_MULTI_UBSPLINE_3D_D_VGL)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,  
   double *val, double *grad, double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vgh,FEVAL_MULTI_UBSPLINE_3D_D_VGH)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad, double *hess);

/////////////////////////////////
// 3D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c,FEVAL_MULTI_UBSPLINE_3D_C)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vg,FEVAL_MULTI_UBSPLINE_3D_C_VG)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z, 
   complex_float *val, complex_float *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vgl,FEVAL_MULTI_UBSPLINE_3D_C_VGL)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vgh,FEVAL_MULTI_UBSPLINE_3D_C_VGH)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *hess);

/////////////////////////////////
// 3D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z,FEVAL_MULTI_UBSPLINE_3D_Z)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vg,FEVAL_MULTI_UBSPLINE_3D_Z_VG)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z, 
   complex_double *val, complex_double *grad);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vgl,FEVAL_MULTI_UBSPLINE_3D_Z_VGL)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *lapl);

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vgh,FEVAL_MULTI_UBSPLINE_3D_Z_VGH)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *hess);


#undef CFUNC
#endif
