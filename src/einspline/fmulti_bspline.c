#include "multi_bspline_create.h"
#include "multi_bspline.h"
#include "fmulti_bspline.h"
#include "config.h"

#ifdef __cplusplus
#define CFUNC "C" /* Avoid name mangling in C++ */
#else
#define CFUNC
#endif


///////////////////////
// Creation routines //
///////////////////////

////////
// 1D //
////////
CFUNC void
F77_FUNC_(fcreate_multi_ubspline_1d_s,FCREATE_MULTI_UBSPLINE_1D_S)
  (double *x0,   double    *x1, int   *num_x, 
   int *x0_code, float *x0_val, int *x1_code, float *x1_val,
   int *num_splines,  multi_UBspline_1d_s **spline)
{
  Ugrid xgrid;
  BCtype_s xBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;

  *spline = create_multi_UBspline_1d_s (xgrid, xBC, *num_splines);
}

CFUNC void
F77_FUNC_(fcreate_multi_ubspline_1d_d,FCREATE_MULTI_UBSPLINE_1D_D)
  (double   *x0, double     *x1, int   *num_x, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *num_splines, multi_UBspline_1d_d **spline)
{
  Ugrid xgrid;
  BCtype_d xBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;

  *spline = create_multi_UBspline_1d_d (xgrid, xBC, *num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_c,FCREATE_MULTI_UBSPLINE_1D_C)
  (double *x0, double *x1, int *num_x, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *num_splines, multi_UBspline_1d_c **spline)
{
  Ugrid xgrid;
  BCtype_c xBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = crealf(*x0_val);
  xBC.lVal_i  = cimagf(*x0_val);
  xBC.rVal_r  = crealf(*x1_val);
  xBC.rVal_i  = cimagf(*x1_val);

  *spline = create_multi_UBspline_1d_c (xgrid, xBC, *num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_1d_z,FCREATE_MULTI_UBSPLINE_1D_Z)
  (double   *x0, double     *x1, int   *num_x, 
   int *x0_code, complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *num_splines, multi_UBspline_1d_z **spline)
{
  Ugrid xgrid;
  BCtype_z xBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = creal(*x0_val);
  xBC.lVal_i  = cimag(*x0_val);
  xBC.rVal_r  = creal(*x1_val);
  xBC.rVal_i  = cimag(*x1_val);

  *spline = create_multi_UBspline_1d_z (xgrid, xBC, *num_splines);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_s,FSET_MULTI_UBSPLINE_1D_S)
  (multi_UBspline_1d_s **spline, int *spline_num, float *data)
{
  set_multi_UBspline_1d_s (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_d,FSET_MULTI_UBSPLINE_1D_D)
  (multi_UBspline_1d_d **spline,  int *spline_num, double *data) 
{
  set_multi_UBspline_1d_d (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_c,FSET_MULTI_UBSPLINE_1D_C)
  (multi_UBspline_1d_c **spline, int *spline_num, complex_float *data)
{
  set_multi_UBspline_1d_c (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_1d_z,FSET_MULTI_UBSPLINE_1D_Z)
  (multi_UBspline_1d_z **spline,  int *spline_num, complex_double *data) 
{
  set_multi_UBspline_1d_z (*spline, *spline_num, data);
}

////////
// 2D //
////////
CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_s,FCREATE_MULTI_UBSPLINE_2D_S)
  (double   *x0, double    *x1, int   *num_x, 
   double   *y0, double    *y1, int   *num_y, 
   int *x0_code, float *x0_val, int *x1_code, float *x1_val,
   int *y0_code, float *y0_val, int *y1_code, float *y1_val,
   int *num_splines, multi_UBspline_2d_s **spline)
{
  Ugrid  xgrid, ygrid;
  BCtype_s xBC, yBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;
  *spline = create_multi_UBspline_2d_s (xgrid, ygrid, xBC, yBC, *num_splines);
}


CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_d,FCREATE_MULTI_UBSPLINE_2D_D)
  (double   *x0, double     *x1, int   *num_x, 
   double   *y0, double     *y1, int   *num_y, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   int *num_splines, multi_UBspline_2d_d **spline)
{
  Ugrid  xgrid, ygrid;
  BCtype_d xBC, yBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal  = *x0_val;
  xBC.rVal  = *x1_val;
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal  = *y0_val;
  yBC.rVal  = *y1_val;
  *spline = create_multi_UBspline_2d_d (xgrid, ygrid, xBC, yBC, *num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_c,FCREATE_MULTI_UBSPLINE_2D_C)
  (double   *x0, double    *x1, int   *num_x, 
   double   *y0, double    *y1, int   *num_y, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, int *y1_code, complex_float *y1_val,
   int *num_splines, multi_UBspline_2d_c **spline)
{
  Ugrid  xgrid, ygrid;
  BCtype_c xBC, yBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
 
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

  *spline = create_multi_UBspline_2d_c (xgrid, ygrid, xBC, yBC, *num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_2d_z,FCREATE_MULTI_UBSPLINE_2D_Z)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   int *x0_code, complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *y0_code, complex_double *y0_val, int *y1_code, complex_double *y1_val,
   int *num_splines, multi_UBspline_2d_z **spline)
{
  Ugrid  xgrid, ygrid;
  BCtype_z xBC, yBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
 
  xBC.lCode = (bc_code) *x0_code;
  xBC.rCode = (bc_code) *x1_code;
  xBC.lVal_r  = crealf(*x0_val);
  xBC.lVal_i  = cimagf(*x0_val);
  xBC.rVal_r  = crealf(*x1_val);
  xBC.rVal_i  = cimagf(*x1_val);
  yBC.lCode = (bc_code) *y0_code;
  yBC.rCode = (bc_code) *y1_code;
  yBC.lVal_r  = creal(*y0_val);
  yBC.lVal_i  = cimag(*y0_val);
  yBC.rVal_r  = creal(*y1_val);
  yBC.rVal_i  = cimag(*y1_val);

  *spline = create_multi_UBspline_2d_z (xgrid, ygrid, xBC, yBC, *num_splines);
}


CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_s,FSET_MULTI_UBSPLINE_2D_S)
  (multi_UBspline_2d_s **spline,  int *spline_num, float *data)
{
  set_multi_UBspline_2d_s (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_d,FSET_MULTI_UBSPLINE_2D_D)
  (multi_UBspline_2d_d **spline, int *spline_num, double *data) 
{
  set_multi_UBspline_2d_d (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_c,FSET_MULTI_UBSPLINE_2D_C)
  (multi_UBspline_2d_c **spline, int *spline_num, complex_float *data)
{
  set_multi_UBspline_2d_c (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_2d_z,FSET_MULTI_UBSPLINE_2D_Z)
  (multi_UBspline_2d_z **spline, int *spline_num, complex_double *data) 
{
  set_multi_UBspline_2d_z (*spline, *spline_num, data);
}


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
   int *num_splines, multi_UBspline_3d_s **spline)
{
  Ugrid  xgrid, ygrid, zgrid;
  BCtype_s xBC, yBC, zBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
  zgrid.start = *z0;
  zgrid.end   = *z1;
  zgrid.num   = *num_z;
 
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
  *spline = create_multi_UBspline_3d_s (xgrid, ygrid, zgrid, xBC, yBC, zBC, 
					*num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_d,FCREATE_MULTI_UBSPLINE_3D_D)
  (double   *x0, double     *x1, int   *num_x, 
   double   *y0, double     *y1, int   *num_y, 
   double   *z0, double     *z1, int   *num_z, 
   int *x0_code, double *x0_val, int *x1_code, double *x1_val,
   int *y0_code, double *y0_val, int *y1_code, double *y1_val,
   int *z0_code, double *z0_val, int *z1_code, double *z1_val,
   int *num_splines, multi_UBspline_3d_d **spline)
{
  Ugrid  xgrid, ygrid, zgrid;
  BCtype_d xBC, yBC, zBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
  zgrid.start = *z0;
  zgrid.end   = *z1;
  zgrid.num   = *num_z;
 
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
  *spline = create_multi_UBspline_3d_d (xgrid, ygrid, zgrid, xBC, yBC, zBC, 
					*num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_c,FCREATE_MULTI_UBSPLINE_3D_C)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   double *z0, double *z1, int *num_z, 
   int *x0_code, complex_float *x0_val, int *x1_code, complex_float *x1_val,
   int *y0_code, complex_float *y0_val, int *y1_code, complex_float *y1_val,
   int *z0_code, complex_float *z0_val, int *z1_code, complex_float *z1_val,
   int *num_splines, multi_UBspline_3d_c **spline)
{
  Ugrid  xgrid, ygrid, zgrid;
  BCtype_c xBC, yBC, zBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
  zgrid.start = *z0;
  zgrid.end   = *z1;
  zgrid.num   = *num_z;
 
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

  *spline = create_multi_UBspline_3d_c (xgrid, ygrid, zgrid, xBC, yBC, zBC, 
					*num_splines);
}

CFUNC void 
F77_FUNC_(fcreate_multi_ubspline_3d_z,FCREATE_MULTI_UBSPLINE_3D_Z)
  (double *x0, double *x1, int *num_x, 
   double *y0, double *y1, int *num_y, 
   double *z0, double *z1, int *num_z, 
   int *x0_code,  complex_double *x0_val, int *x1_code, complex_double *x1_val,
   int *y0_code,  complex_double *y0_val, int *y1_code, complex_double *y1_val,
   int *z0_code,  complex_double *z0_val, int *z1_code, complex_double *z1_val,
   int *num_splines, multi_UBspline_3d_z **spline)
{
  Ugrid  xgrid, ygrid, zgrid;
  BCtype_z xBC, yBC, zBC;
  xgrid.start = *x0;
  xgrid.end   = *x1;
  xgrid.num   = *num_x;
  ygrid.start = *y0;
  ygrid.end   = *y1;
  ygrid.num   = *num_y;
  zgrid.start = *z0;
  zgrid.end   = *z1;
  zgrid.num   = *num_z;
 
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

  *spline = create_multi_UBspline_3d_z (xgrid, ygrid, zgrid, xBC, yBC, zBC, 
					*num_splines);
}


CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_s,FSET_MULTI_UBSPLINE_3D_S)
  (multi_UBspline_3d_s **spline, int *spline_num, float *data)
{
  set_multi_UBspline_3d_s (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_d,FSET_MULTI_UBSPLINE_3D_D)
  (multi_UBspline_3d_d **spline, int *spline_num, double *data) 
{
  set_multi_UBspline_3d_d (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_c,FSET_MULTI_UBSPLINE_3D_C)
  (multi_UBspline_3d_c **spline, int *spline_num, complex_float *data)
{
  set_multi_UBspline_3d_c (*spline, *spline_num, data);
}

CFUNC void 
F77_FUNC_(fset_multi_ubspline_3d_z,FSET_MULTI_UBSPLINE_3D_Z)
  (multi_UBspline_3d_z **spline, int *spline_num, complex_double *data) 
{
  set_multi_UBspline_3d_z (*spline, *spline_num, data);
}


/////////////////////////
// Evaluation routines //
/////////////////////////

//////////////////////////////
// 1D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_s,FEVAL_MULTI_UBSPLINE_1D_S)
  (multi_UBspline_1d_s **spline, double *x, float *val)
{
  eval_multi_UBspline_1d_s (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_s_vg,FEVAL_MULTI_UBSPLINE_1D_S_VG)
  (multi_UBspline_1d_s **spline, double *x, float *val, float *grad)
{
  eval_multi_UBspline_1d_s_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_s_vgl,FEVAL_MULTI_UBSPLINE_1D_S_VGL)
  (multi_UBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *lapl)
{
  eval_multi_UBspline_1d_s_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_s_vgh,FEVAL_MULTI_UBSPLINE_1D_S_VGH)
  (multi_UBspline_1d_s **spline, double *x, 
   float *val, float *grad, float *hess)
{
  eval_multi_UBspline_1d_s_vgh (*spline, *x, val, grad, hess);
}

//////////////////////////////
// 1D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d,FEVAL_MULTI_UBSPLINE_1D_D)
  (multi_UBspline_1d_d **spline, double *x, double *val)
{
  eval_multi_UBspline_1d_d (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vg,FEVAL_MULTI_UBSPLINE_1D_D_VG)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad)
{
  eval_multi_UBspline_1d_d_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vgl,FEVAL_MULTI_UBSPLINE_1D_D_VGL)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *lapl)
{
  eval_multi_UBspline_1d_d_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_d_vgh,FEVAL_MULTI_UBSPLINE_1D_D_VGH)
  (multi_UBspline_1d_d **spline, double *x, 
   double *val, double *grad, double *hess)
{
  eval_multi_UBspline_1d_d_vgh (*spline, *x, val, grad, hess);
}

/////////////////////////////////
// 1D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c,FEVAL_MULTI_UBSPLINE_1D_C)
  (multi_UBspline_1d_c **spline, double *x, complex_float *val)
{
  eval_multi_UBspline_1d_c (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vg,FEVAL_MULTI_UBSPLINE_1D_C_VG)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad)
{
  eval_multi_UBspline_1d_c_vg (*spline, *x, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vgl,FEVAL_MULTI_UBSPLINE_1D_C_VGL)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_multi_UBspline_1d_c_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_c_vgh,FEVAL_MULTI_UBSPLINE_1D_C_VGH)
  (multi_UBspline_1d_c **spline, double *x, 
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_multi_UBspline_1d_c_vgh (*spline, *x, val, grad, hess);
}

/////////////////////////////////
// 1D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z,FEVAL_MULTI_UBSPLINE_1D_Z)
  (multi_UBspline_1d_z **spline, double *x, complex_double *val)
{
  eval_multi_UBspline_1d_z (*spline, *x, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vg,FEVAL_MULTI_UBSPLINE_1D_Z_VG)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad)
{
  eval_multi_UBspline_1d_z_vg (*spline, *x, val, grad);
}
 
CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vgl,FEVAL_MULTI_UBSPLINE_1D_Z_VGL)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_multi_UBspline_1d_z_vgl (*spline, *x, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_1d_z_vgh,FEVAL_MULTI_UBSPLINE_1D_Z_VGH)
  (multi_UBspline_1d_z **spline, double *x, 
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_multi_UBspline_1d_z_vgh (*spline, *x, val, grad, hess);
}

//////////////////////////////
// 2D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s,FEVAL_MULTI_UBSPLINE_2D_S)
  (multi_UBspline_2d_s **spline, double *x, double *y, float *val)
{
  eval_multi_UBspline_2d_s (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vg,FEVAL_MULTI_UBSPLINE_2D_S_VG)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad)
{
  eval_multi_UBspline_2d_s_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vgl,FEVAL_MULTI_UBSPLINE_2D_S_VGL)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float* lapl)
{
  eval_multi_UBspline_2d_s_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_s_vgh,FEVAL_MULTI_UBSPLINE_2D_S_VGH)
  (multi_UBspline_2d_s **spline, double *x, double *y, 
   float *val, float *grad, float *hess)
{
  eval_multi_UBspline_2d_s_vgh (*spline, *x, *y, val, grad, hess);
}

//////////////////////////////
// 2D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d,FEVAL_MULTI_UBSPLINE_2D_D)
  (multi_UBspline_2d_d **spline, double *x, double *y, double *val)
{
  eval_multi_UBspline_2d_d (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vg,FEVAL_MULTI_UBSPLINE_2D_D_VG)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad)
{
  eval_multi_UBspline_2d_d_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vgl,FEVAL_MULTI_UBSPLINE_2D_D_VGL)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *lapl)
{
  eval_multi_UBspline_2d_d_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_d_vgh,FEVAL_MULTI_UBSPLINE_2D_D_VGH)
  (multi_UBspline_2d_d **spline, double *x, double *y, 
   double *val, double *grad, double *hess)
{
  eval_multi_UBspline_2d_d_vgl (*spline, *x, *y, val, grad, hess);
}

/////////////////////////////////
// 2D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c,FEVAL_MULTI_UBSPLINE_2D_C)
  (multi_UBspline_2d_c **spline, double *x, double *y, complex_float *val)
{
  eval_multi_UBspline_2d_c (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vg,FEVAL_MULTI_UBSPLINE_2D_C_VG)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad)
{
  eval_multi_UBspline_2d_c_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vgl,FEVAL_MULTI_UBSPLINE_2D_C_VGL)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_multi_UBspline_2d_c_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_c_vgh,FEVAL_MULTI_UBSPLINE_2D_C_VGH)
  (multi_UBspline_2d_c **spline, double *x, double *y, 
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_multi_UBspline_2d_c_vgh (*spline, *x, *y, val, grad, hess);
}

/////////////////////////////////
// 2D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z,FEVAL_MULTI_UBSPLINE_2D_Z)
  (multi_UBspline_2d_z **spline, double *x, double *y, complex_double *val)
{
  eval_multi_UBspline_2d_z (*spline, *x, *y, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vg,FEVAL_MULTI_UBSPLINE_2D_Z_VG)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad)
{
  eval_multi_UBspline_2d_z_vg (*spline, *x, *y, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vgl,FEVAL_MULTI_UBSPLINE_2D_Z_VGL)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_multi_UBspline_2d_z_vgl (*spline, *x, *y, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_2d_z_vgh,FEVAL_MULTI_UBSPLINE_2D_Z_VGH)
  (multi_UBspline_2d_z **spline, double *x, double *y, 
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_multi_UBspline_2d_z_vgh (*spline, *x, *y, val, grad, hess);
}



//////////////////////////////
// 3D single-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s,FEVAL_MULTI_UBSPLINE_3D_S)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val)
{
  eval_multi_UBspline_3d_s (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vg,FEVAL_MULTI_UBSPLINE_3D_S_VG)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad)
{
  eval_multi_UBspline_3d_s_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vgl,FEVAL_MULTI_UBSPLINE_3D_S_VGL)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z,
   float *val, float *grad, float* lapl)
{
  eval_multi_UBspline_3d_s_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_s_vgh,FEVAL_MULTI_UBSPLINE_3D_S_VGH)
  (multi_UBspline_3d_s **spline, double *x, double *y, double *z, 
   float *val, float *grad, float *hess)
{
  eval_multi_UBspline_3d_s_vgh (*spline, *x, *y, *z, val, grad, hess);
}

//////////////////////////////
// 3D double-precision real //
//////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d,FEVAL_MULTI_UBSPLINE_3D_D)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val)
{
  eval_multi_UBspline_3d_d (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vg,FEVAL_MULTI_UBSPLINE_3D_D_VG)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad)
{
  eval_multi_UBspline_3d_d_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vgl,FEVAL_MULTI_UBSPLINE_3D_D_VGL)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,  
   double *val, double *grad, double *lapl)
{
  eval_multi_UBspline_3d_d_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_d_vgh,FEVAL_MULTI_UBSPLINE_3D_D_VGH)
  (multi_UBspline_3d_d **spline, double *x, double *y, double *z,
   double *val, double *grad, double *hess)
{
  eval_multi_UBspline_3d_d_vgh (*spline, *x, *y, *z, val, grad, hess);
}

/////////////////////////////////
// 3D single-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c,FEVAL_MULTI_UBSPLINE_3D_C)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val)
{
  eval_multi_UBspline_3d_c (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vg,FEVAL_MULTI_UBSPLINE_3D_C_VG)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z, 
   complex_float *val, complex_float *grad)
{
  eval_multi_UBspline_3d_c_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vgl,FEVAL_MULTI_UBSPLINE_3D_C_VGL)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *lapl)
{
  eval_multi_UBspline_3d_c_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_c_vgh,FEVAL_MULTI_UBSPLINE_3D_C_VGH)
  (multi_UBspline_3d_c **spline, double *x, double *y, double *z,
   complex_float *val, complex_float *grad, complex_float *hess)
{
  eval_multi_UBspline_3d_c_vgh (*spline, *x, *y, *z, val, grad, hess);
}

/////////////////////////////////
// 3D double-precision complex //
/////////////////////////////////
CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z,FEVAL_MULTI_UBSPLINE_3D_Z)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val)
{
  eval_multi_UBspline_3d_z (*spline, *x, *y, *z, val);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vg,FEVAL_MULTI_UBSPLINE_3D_Z_VG)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z, 
   complex_double *val, complex_double *grad)
{
  eval_multi_UBspline_3d_z_vg (*spline, *x, *y, *z, val, grad);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vgl,FEVAL_MULTI_UBSPLINE_3D_Z_VGL)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *lapl)
{
  eval_multi_UBspline_3d_z_vgl (*spline, *x, *y, *z, val, grad, lapl);
}

CFUNC void
F77_FUNC_(feval_multi_ubspline_3d_z_vgh,FEVAL_MULTI_UBSPLINE_3D_Z_VGH)
  (multi_UBspline_3d_z **spline, double *x, double *y, double *z,
   complex_double *val, complex_double *grad, complex_double *hess)
{
  eval_multi_UBspline_3d_z_vgh (*spline, *x, *y, *z, val, grad, hess);
}


