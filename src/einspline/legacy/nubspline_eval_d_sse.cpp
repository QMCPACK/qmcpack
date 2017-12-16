/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//                                                                         //
//  This program is free software; you can redistribute it and/or modify   //
//  it under the terms of the GNU General Public License as published by   //
//  the Free Software Foundation; either version 2 of the License, or      //
//  (at your option) any later version.                                    //
//                                                                         //
//  This program is distributed in the hope that it will be useful,        //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//  GNU General Public License for more details.                           //
//                                                                         //
//  You should have received a copy of the GNU General Public License      //
//  along with this program; if not, write to the Free Software            //
//  Foundation, Inc., 51 Franklin Street, Fifth Floor,                     //
//  Boston, MA  02110-1301  USA                                            //
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <config.h>
#include "bspline_base.h"
#include "nubspline_structs.h"
#include "nubspline_eval_d.h"

#ifdef HAVE_SSE2
#include <xmmintrin.h>
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif

#ifdef HAVE_SSE3
#define _MM_DDOT4_PD(a0, a1, a2, a3, b0, b1, b2, b3, r)               \
do {                                                                  \
   __m128d t0 = _mm_add_pd(_mm_mul_pd (a0, b0),_mm_mul_pd (a1, b1));  \
   __m128d t1 = _mm_add_pd(_mm_mul_pd (a2, b2),_mm_mul_pd (a3, b3));  \
   r = _mm_hadd_pd (t0, t1);                                          \
 } while(0);
#define _MM_DOT4_PD(a0, a1, b0, b1, p)                                \
do {                                                                  \
  __m128d t0 = _mm_add_pd(_mm_mul_pd (a0, b0),_mm_mul_pd (a1, b1));   \
  __m128d t1 = _mm_hadd_pd (t0,t0);                                   \
  _mm_store_sd (&(p), t1);                                            \
 } while (0);
#else
#define _MM_DDOT4_PD(a0, a1, a2, a3, b0, b1, b2, b3, r)               \
do {                                                                  \
   __m128d t0 = _mm_add_pd(_mm_mul_pd (a0, b0),_mm_mul_pd (a1, b1));  \
   __m128d t1 = _mm_add_pd(_mm_mul_pd (a2, b2),_mm_mul_pd (a3, b3));  \
   r = _mm_add_pd(_mm_unpacklo_pd(t0,t1),_mm_unpackhi_pd(t0,t1));     \
 } while(0);
#define _MM_DOT4_PD(a0, a1, b0, b1, p)                                \
do {                                                                  \
  __m128d t0 = _mm_add_pd(_mm_mul_pd (a0, b0),_mm_mul_pd (a1, b1));   \
  __m128d t1 =                                                        \
      _mm_add_pd (_mm_unpacklo_pd(t0,t0), _mm_unpackhi_pd(t0,t0));    \
  _mm_store_sd (&(p), t1);                                            \
 } while (0);
#endif


/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void
eval_NUBspline_1d_d (NUBspline_1d_d * restrict spline,
                     double x, double* restrict val)
{
  double bfuncs[4];
  int i = get_NUBasis_funcs_d (spline->x_basis, x, bfuncs);
  double* restrict coefs = spline->coefs;
  *val = (coefs[i+0]*bfuncs[0] +coefs[i+1]*bfuncs[1] +
          coefs[i+2]*bfuncs[2] +coefs[i+3]*bfuncs[3]);
}

/* Value and first derivative */
void
eval_NUBspline_1d_d_vg (NUBspline_1d_d * restrict spline, double x,
                        double* restrict val, double* restrict grad)
{
  double bfuncs[4], dbfuncs[4];
  int i = get_NUBasis_dfuncs_d (spline->x_basis, x, bfuncs, dbfuncs);
  double* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]* bfuncs[0] + coefs[i+1]* bfuncs[1] +
           coefs[i+2]* bfuncs[2] + coefs[i+3]* bfuncs[3]);
  *grad = (coefs[i+0]*dbfuncs[0] + coefs[i+1]*dbfuncs[1] +
           coefs[i+2]*dbfuncs[2] + coefs[i+3]*dbfuncs[3]);
}

/* Value, first derivative, and second derivative */
void
eval_NUBspline_1d_d_vgl (NUBspline_1d_d * restrict spline, double x,
                         double* restrict val, double* restrict grad,
                         double* restrict lapl)
{
  double bfuncs[4], dbfuncs[4], d2bfuncs[4];
  int i = get_NUBasis_d2funcs_d (spline->x_basis, x, bfuncs, dbfuncs, d2bfuncs);
  double* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]*  bfuncs[0] + coefs[i+1]*  bfuncs[1] +
           coefs[i+2]*  bfuncs[2] + coefs[i+3]*  bfuncs[3]);
  *grad = (coefs[i+0]* dbfuncs[0] + coefs[i+1]* dbfuncs[1] +
           coefs[i+2]* dbfuncs[2] + coefs[i+3]* dbfuncs[3]);
  *lapl = (coefs[i+0]*d2bfuncs[0] + coefs[i+1]*d2bfuncs[1] +
           coefs[i+2]*d2bfuncs[2] + coefs[i+3]*d2bfuncs[3]);
}

void
eval_NUBspline_1d_d_vgh (NUBspline_1d_d * restrict spline, double x,
                         double* restrict val, double* restrict grad,
                         double* restrict hess)
{
  eval_NUBspline_1d_d_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void
eval_NUBspline_2d_d (NUBspline_2d_d * restrict spline,
                     double x, double y, double* restrict val)
{
  __m128d a01, b01, bP01,
          a23, b23, bP23,
          tmp0, tmp1, tmp2, tmp3;
  int ix = get_NUBasis_funcs_sse_d (spline->x_basis, x, &a01, &a23);
  int iy = get_NUBasis_funcs_sse_d (spline->y_basis, y, &b01, &b23);
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0));
  tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0));
  tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  tmp0 = _mm_loadu_pd (P(2,0));
  tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0));
  tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
#undef P
}


/* Value and gradient */
void
eval_NUBspline_2d_d_vg (NUBspline_2d_d * restrict spline,
                        double x, double y,
                        double* restrict val, double* restrict grad)
{
  __m128d a01, b01, da01, db01, bP01, dbP01,
          a23, b23, da23, db23, bP23, dbP23,
          tmp0, tmp1, tmp2, tmp3;
  int ix = get_NUBasis_dfuncs_sse_d (spline->x_basis, x,
                                     &a01, &a23, &da01, &da23);
  int iy = get_NUBasis_dfuncs_sse_d (spline->y_basis, y,
                                     &b01, &b23, &db01, &db23);
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0));
  tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0));
  tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  tmp0 = _mm_loadu_pd (P(2,0));
  tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0));
  tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
#undef P
}

/* Value, gradient, and laplacian */
void
eval_NUBspline_2d_d_vgl (NUBspline_2d_d * restrict spline,
                         double x, double y, double* restrict val,
                         double* restrict grad, double* restrict lapl)
{
  __m128d a01, b01, da01, db01, d2a01, d2b01,
          a23, b23, da23, db23, d2a23, d2b23,
          bP01, dbP01, d2bP01,
          bP23, dbP23, d2bP23,
          tmp0, tmp1, tmp2, tmp3;
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x,
                                      &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y,
                                      &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0));
  tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0));
  tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP01);
  tmp0 = _mm_loadu_pd (P(2,0));
  tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0));
  tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
  // Compute laplacian
  double d2x, d2y;
  _MM_DOT4_PD (d2a01, d2a23,   bP01,   bP23, d2x);
  _MM_DOT4_PD (  a01,   a23, d2bP01, d2bP23, d2y);
  *lapl = d2x + d2y;
#undef P
}

/* Value, gradient, and Hessian */
void
eval_NUBspline_2d_d_vgh (NUBspline_2d_d * restrict spline,
                         double x, double y, double* restrict val,
                         double* restrict grad, double* restrict hess)
{
  __m128d a01, b01, da01, db01, d2a01, d2b01,
          a23, b23, da23, db23, d2a23, d2b23,
          bP01, dbP01, d2bP01,
          bP23, dbP23, d2bP23,
          tmp0, tmp1, tmp2, tmp3;
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x,
                                      &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y,
                                      &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0));
  tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0));
  tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP01);
  tmp0 = _mm_loadu_pd (P(2,0));
  tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0));
  tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
  // Compute hessian
  _MM_DOT4_PD (d2a01, d2a23,   bP01,   bP23, hess[0]);
  _MM_DOT4_PD (  a01,   a23, d2bP01, d2bP23, hess[3]);
  _MM_DOT4_PD ( da01,  da23,  dbP01,  dbP23, hess[1]);
  hess[2] = hess[1];
#undef P
}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void
eval_NUBspline_3d_d (NUBspline_3d_d * restrict spline,
                     double x, double y, double z,
                     double* restrict val)
{
  __m128d
  a01, b01, c01, a23, b23, c23,
       cP[8], dcP[8], d2cP[8],
       bcP01, dbcP01, bdcP01,
       bcP23, dbcP23, bdcP23,
       tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int ix = get_NUBasis_funcs_sse_d (spline->x_basis, x, &a01, &a23);
  int iy = get_NUBasis_funcs_sse_d (spline->y_basis, y, &b01, &b23);
  int iz = get_NUBasis_funcs_sse_d (spline->z_basis, z, &c01, &c23);
  int xs   = spline->x_stride;
  int ys   = spline->y_stride;
  int ysb  = ys+2;
  int ys2  = 2*ys;
  int ys2b = 2*ys + 2;
  int ys3  = 3*ys;
  int ys3b = 3*ys + 2;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  double *p = P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  p = P(0,0,0);
  // 1st eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  // 2nd eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  p += xs;
  // 3rd eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  // 4th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  p += xs;
  // 5th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  // 6th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  p += xs;
  // 7th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  // 8th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  // Now compute bcP, dbcP, bdcP and products
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[0], cP[1], cP[2], cP[3], bcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[4], cP[5], cP[6], cP[7], bcP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);
#undef P
}

/* Value and gradient */
void
eval_NUBspline_3d_d_vg (NUBspline_3d_d * restrict spline,
                        double x, double y, double z,
                        double* restrict val, double* restrict grad)
{
  __m128d
  a01, b01, c01, da01, db01, dc01,
       a23, b23, c23, da23, db23, dc23,
       cP[8], dcP[8], d2cP[8],
       bcP01, dbcP01, bdcP01,
       bcP23, dbcP23, bdcP23,
       tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int ix = get_NUBasis_dfuncs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23);
  int iy = get_NUBasis_dfuncs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23);
  int iz = get_NUBasis_dfuncs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23);
  int xs   = spline->x_stride;
  int ys   = spline->y_stride;
  int ysb  = ys+2;
  int ys2  = 2*ys;
  int ys2b = 2*ys + 2;
  int ys3  = 3*ys;
  int ys3b = 3*ys + 2;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  double *p = P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  p = P(0,0,0);
  // 1st eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);
  // 2nd eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);
  p += xs;
  // 3rd eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);
  // 4th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);
  p += xs;
  // 5th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);
  // 6th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);
  p += xs;
  // 7th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);
  // 8th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[7]);
  // Now compute bcP, dbcP, bdcP and products
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[0], cP[1], cP[2], cP[3], bcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[4], cP[5], cP[6], cP[7], bcP23);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[0], cP[1], cP[2], cP[3], dbcP01);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[4], cP[5], cP[6], cP[7], dbcP23);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[0], dcP[1], dcP[2], dcP[3], bdcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[4], dcP[5], dcP[6], dcP[7], bdcP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bcP01, bcP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbcP01, dbcP23, grad[1]);
  _MM_DOT4_PD (a01, a23, bdcP01, bdcP23, grad[2]);
#undef P
}



/* Value, gradient, and laplacian */
void
eval_NUBspline_3d_d_vgl (NUBspline_3d_d * restrict spline,
                         double x, double y, double z,
                         double* restrict val, double* restrict grad, double* restrict lapl)
{
  __m128d
  a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
       a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
       cP[8], dcP[8], d2cP[8],
       bcP01, dbcP01, bdcP01, d2bcP01, dbdcP01, bd2cP01,
       bcP23, dbcP23, bdcP23, d2bcP23, dbdcP23, bd2cP23,
       tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int iz = get_NUBasis_d2funcs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23, &d2c01, &d2c23);
  int xs   = spline->x_stride;
  int ys   = spline->y_stride;
  int ysb  = ys+2;
  int ys2  = 2*ys;
  int ys2b = 2*ys + 2;
  int ys3  = 3*ys;
  int ys3b = 3*ys + 2;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  double *p = P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  p = P(0,0,0);
  // 1st eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[0]);
  // 2nd eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[1]);
  p += xs;
  // 3rd eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[2]);
  // 4th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[3]);
  p += xs;
  // 5th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[4]);
  // 6th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[5]);
  p += xs;
  // 7th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[6]);
  // 8th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[7]);
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[0], cP[1], cP[2], cP[3], bcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[4], cP[5], cP[6], cP[7], bcP23);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[0], cP[1], cP[2], cP[3], dbcP01);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[4], cP[5], cP[6], cP[7], dbcP23);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[0], dcP[1], dcP[2], dcP[3], bdcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[4], dcP[5], dcP[6], dcP[7], bdcP23);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cP[0], cP[1], cP[2], cP[3], d2bcP01);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cP[4], cP[5], cP[6], cP[7], d2bcP23);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cP[0], d2cP[1], d2cP[2], d2cP[3], bd2cP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cP[4], d2cP[5], d2cP[6], d2cP[7], bd2cP23);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcP[0], dcP[1], dcP[2], dcP[3], dbdcP01);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcP[4], dcP[5], dcP[6], dcP[7], dbdcP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bcP01, bcP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbcP01, dbcP23, grad[1]);
  _MM_DOT4_PD (a01, a23, bdcP01, bdcP23, grad[2]);
  // Compute Laplacian
  double d2x, d2y, d2z;
  _MM_DOT4_PD (d2a01, d2a23, bcP01, bcP23, d2x);
  _MM_DOT4_PD (a01, a23, d2bcP01, d2bcP23, d2y);
  _MM_DOT4_PD (a01, a23, bd2cP01, bd2cP23, d2z);
  *lapl = d2x + d2y + d2z;
#undef P
}





/* Value, gradient, and Hessian */
void
eval_NUBspline_3d_d_vgh (NUBspline_3d_d * restrict spline,
                         double x, double y, double z,
                         double* restrict val, double* restrict grad, double* restrict hess)
{
  __m128d
  a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
       a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
       cP[8], dcP[8], d2cP[8],
       bcP01, dbcP01, bdcP01, d2bcP01, dbdcP01, bd2cP01,
       bcP23, dbcP23, bdcP23, d2bcP23, dbdcP23, bd2cP23,
       tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int iz = get_NUBasis_d2funcs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23, &d2c01, &d2c23);
  int xs   = spline->x_stride;
  int ys   = spline->y_stride;
  int ysb  = ys+2;
  int ys2  = 2*ys;
  int ys2b = 2*ys + 2;
  int ys3  = 3*ys;
  int ys3b = 3*ys + 2;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  double *p = P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p     ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+2   ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys  ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ysb ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2b), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3 ), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3b), _MM_HINT_T0);
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  p = P(0,0,0);
  // 1st eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[0]);
  // 2nd eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[1]);
  p += xs;
  // 3rd eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[2]);
  // 4th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[3]);
  p += xs;
  // 5th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[4]);
  // 6th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[5]);
  p += xs;
  // 7th eighth
  tmp0    = _mm_loadu_pd (p    );
  tmp1    = _mm_loadu_pd (p+2  );
  tmp2    = _mm_loadu_pd (p+ys );
  tmp3    = _mm_loadu_pd (p+ysb);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[6]);
  // 8th eighth
  tmp0    = _mm_loadu_pd (p+ys2 );
  tmp1    = _mm_loadu_pd (p+ys2b);
  tmp2    = _mm_loadu_pd (p+ys3 );
  tmp3    = _mm_loadu_pd (p+ys3b);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[7]);
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[0], cP[1], cP[2], cP[3], bcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, cP[4], cP[5], cP[6], cP[7], bcP23);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[0], cP[1], cP[2], cP[3], dbcP01);
  _MM_DDOT4_PD (db01, db23, db01, db23, cP[4], cP[5], cP[6], cP[7], dbcP23);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[0], dcP[1], dcP[2], dcP[3], bdcP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcP[4], dcP[5], dcP[6], dcP[7], bdcP23);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cP[0], cP[1], cP[2], cP[3], d2bcP01);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cP[4], cP[5], cP[6], cP[7], d2bcP23);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cP[0], d2cP[1], d2cP[2], d2cP[3], bd2cP01);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cP[4], d2cP[5], d2cP[6], d2cP[7], bd2cP23);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcP[0], dcP[1], dcP[2], dcP[3], dbdcP01);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcP[4], dcP[5], dcP[6], dcP[7], dbdcP23);
  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bcP01, bcP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbcP01, dbcP23, grad[1]);
  _MM_DOT4_PD (a01, a23, bdcP01, bdcP23, grad[2]);
  // Compute hessian
  // d2x
  _MM_DOT4_PD (d2a01, d2a23, bcP01, bcP23, hess[0]);
  // d2y
  _MM_DOT4_PD (a01, a23, d2bcP01, d2bcP23, hess[4]);
  // d2z
  _MM_DOT4_PD (a01, a23, bd2cP01, bd2cP23, hess[8]);
  // dx dy
  _MM_DOT4_PD (da01, da23, dbcP01, dbcP23, hess[1]);
  // dx dz
  _MM_DOT4_PD (da01, da23, bdcP01, bdcP23, hess[2]);
  // dy dz
  _MM_DOT4_PD (a01, a23, dbdcP01, dbdcP23, hess[5]);
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];
#undef P
}

#undef _MM_DDOT4_PD
#undef _MM_DOT4_PD
