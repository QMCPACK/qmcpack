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

#ifndef NUBSPLINE_EVAL_SSE_C_H
#define NUBSPLINE_EVAL_SSE_C_H

#include <math.h>
#include <stdio.h>
#include "nubspline_structs.h"

#ifdef HAVE_SSE
#include <xmmintrin.h>
#endif

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

/// SSE3 adds "horizontal add" instructions, which makes things
/// simpler and faster
#ifdef HAVE_SSE3

#include <pmmintrin.h>

#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 _r0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
  __m128 _r1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
  r = _mm_hadd_ps (_r0, _r1);                                         \
 } while (0);
#define _MM_DOT4_PS(_A, _B, _p)                                       \
do {                                                                \
  __m128 t  = _mm_mul_ps (_A, _B);                                    \
  __m128 t1 = _mm_hadd_ps (t,t);                                    \
  __m128 r  = _mm_hadd_ps (t1, t1);                                 \
  _mm_store_ss (&(_p), r);                                           \
} while(0);
#else
// Use plain-old SSE instructions
#define _MM_MATVEC4_PS(_M0, _M1, _M2, _M3, _v, _r)                        \
do {                                                                \
  __m128 _r0 = _mm_mul_ps (_M0, _v);                                   \
  __m128 _r1 = _mm_mul_ps (_M1, _v);				    \
  __m128 _r2 = _mm_mul_ps (_M2, _v);                                   \
  __m128 _r3 = _mm_mul_ps (_M3, _v);				    \
  _MM_TRANSPOSE4_PS (_r0, _r1, _r2, _r3);                               \
  _r = _mm_add_ps (_mm_add_ps (r0, r1), _mm_add_ps (r2, r3));        \
 } while (0);
#define _MM_DOT4_PS(_A, _B, _p)                                        \
do {                                                                \
  __m128 t    = _mm_mul_ps (_A, _B);                                  \
  __m128 alo  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(0,1,0,1));	    \
  __m128 ahi  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(2,3,2,3));	    \
  __m128 a    = _mm_add_ps (alo, ahi);                              \
  __m128 rlo  = _mm_shuffle_ps (a, a, _MM_SHUFFLE(0,0,0,0));	     \
  __m128 rhi  = _mm_shuffle_ps (a, a, _MM_SHUFFLE(1,1,1,1));	     \
  __m128 r    = _mm_add_ps (rlo, rhi);                              \
  _mm_store_ss (&(_p), r);                                           \
} while(0);
#endif

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_NUBspline_1d_c (NUBspline_1d_c * restrict spline, 
		     double x, complex_float* restrict val)
{
  float bfuncs[4];
  int i = get_NUBasis_funcs_s (spline->x_basis, x, bfuncs);
  complex_float* restrict coefs = spline->coefs;
  *val = (coefs[i+0]*bfuncs[0] +coefs[i+1]*bfuncs[1] +
	  coefs[i+2]*bfuncs[2] +coefs[i+3]*bfuncs[3]);
}

/* Value and first derivative */
inline void
eval_NUBspline_1d_c_vg (NUBspline_1d_c * restrict spline, double x, 
			complex_float* restrict val, complex_float* restrict grad)
{
  float bfuncs[4], dbfuncs[4];
  int i = get_NUBasis_dfuncs_s (spline->x_basis, x, bfuncs, dbfuncs);
  complex_float* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]* bfuncs[0] + coefs[i+1]* bfuncs[1] +
	   coefs[i+2]* bfuncs[2] + coefs[i+3]* bfuncs[3]);
  *grad = (coefs[i+0]*dbfuncs[0] + coefs[i+1]*dbfuncs[1] +
	   coefs[i+2]*dbfuncs[2] + coefs[i+3]*dbfuncs[3]);
}

/* Value, first derivative, and second derivative */
inline void
eval_NUBspline_1d_c_vgl (NUBspline_1d_c * restrict spline, double x, 
			complex_float* restrict val, complex_float* restrict grad,
			complex_float* restrict lapl)
{
  float bfuncs[4], dbfuncs[4], d2bfuncs[4];
  int i = get_NUBasis_d2funcs_s (spline->x_basis, x, bfuncs, dbfuncs, d2bfuncs);
  complex_float* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]*  bfuncs[0] + coefs[i+1]*  bfuncs[1] +
	   coefs[i+2]*  bfuncs[2] + coefs[i+3]*  bfuncs[3]);
  *grad = (coefs[i+0]* dbfuncs[0] + coefs[i+1]* dbfuncs[1] +
	   coefs[i+2]* dbfuncs[2] + coefs[i+3]* dbfuncs[3]);
  *lapl = (coefs[i+0]*d2bfuncs[0] + coefs[i+1]*d2bfuncs[1] +
	   coefs[i+2]*d2bfuncs[2] + coefs[i+3]*d2bfuncs[3]);

}

inline void
eval_NUBspline_1d_c_vgh (NUBspline_1d_c * restrict spline, double x, 
			complex_float* restrict val, complex_float* restrict grad,
			complex_float* restrict hess)
{
  eval_NUBspline_1d_c_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_NUBspline_2d_c (NUBspline_2d_c * restrict spline, 
		     double x, double y, complex_float* restrict val)
{
  __m128 a, b, bPr, bPi, r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;
  complex_float* restrict coefs = spline->coefs;
  int ix = get_NUBasis_funcs_sse_s (spline->x_basis, x, &a);
  int iy = get_NUBasis_funcs_sse_s (spline->y_basis, y, &b);

  int xs = spline->x_stride;
  int xs2 = 2*xs;

  #define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  float* restrict p = (float*)P(0,0);
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);

  p = (float *)P(0,0);
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));

  _MM_MATVEC4_PS (r0, r1, r2, r3,   b,   bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   b,   bPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;

  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);

#undef P
}


/* Value and gradient */
inline void
eval_NUBspline_2d_c_vg (NUBspline_2d_c * restrict spline, 
		       double x, double y, 
		       complex_float* restrict val, complex_float* restrict grad)
{
  __m128 a, b, da, db, 
    bPr, dbPr, bPi, dbPi, r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;
  int ix = get_NUBasis_dfuncs_sse_s (spline->x_basis, x, &a, &da);
  int iy = get_NUBasis_dfuncs_sse_s (spline->y_basis, y, &b, &db);
  
  complex_float* restrict coefs = spline->coefs;

  int xs = spline->x_stride;
  int xs2 = 2*xs;

  #define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  float* restrict p = (float*)P(0,0);
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);

  p = (float *)P(0,0);
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));

  _MM_MATVEC4_PS (r0, r1, r2, r3,   b,   bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   b,   bPi);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  db,  dbPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  db,  dbPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;

  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bPr, *gradr0);
  _MM_DOT4_PS (da, bPi, *gradr0);
  _MM_DOT4_PS (a, dbPr, *gradi1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_NUBspline_2d_c_vgl (NUBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, complex_float* restrict lapl)
{
  __m128 a, b, da, db, d2a, d2b, 
    bPr, dbPr, d2bPr, bPi, dbPi, d2bPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;
  int ix = get_NUBasis_d2funcs_sse_s (spline->x_basis, x, &a, &da, &d2a);
  int iy = get_NUBasis_d2funcs_sse_s (spline->y_basis, y, &b, &db, &d2b);

  complex_float* restrict coefs = spline->coefs;

  int xs = spline->x_stride;
  int xs2 = 2*xs;

  #define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  float* restrict p = (float*)P(0,0);
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);

  p = (float *)P(0,0);
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));

  _MM_MATVEC4_PS (r0, r1, r2, r3,   b,   bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   b,   bPi);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  db,  dbPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  db,  dbPi);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2b, d2bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2b, d2bPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;

  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bPr, *gradr0);
  _MM_DOT4_PS (da, bPi, *gradr0);
  _MM_DOT4_PS (a, dbPr, *gradi1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
  // Compute laplacian
  float d2x_r, d2y_r, d2x_i, d2y_i;
  _MM_DOT4_PS (d2a, bPr, d2x_r);
  _MM_DOT4_PS (d2a, bPi, d2x_i);
  _MM_DOT4_PS (a, d2bPr, d2y_r);
  _MM_DOT4_PS (a, d2bPi, d2y_i);
#ifdef __cplusplus
  *lapl = std::complex<float>(d2x_r + d2y_r,
			      d2x_i + d2y_i);
#else
  *lapl = (d2x_r + d2y_r) + 1.0if*(d2x_i + d2y_i);
#endif
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_NUBspline_2d_c_vgh (NUBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, complex_float* restrict hess)
{
  __m128 a, b, da, db, d2a, d2b, 
    bPr, dbPr, d2bPr, bPi, dbPi, d2bPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;

  int ix = get_NUBasis_d2funcs_sse_s (spline->x_basis, x, &a, &da, &d2a);
  int iy = get_NUBasis_d2funcs_sse_s (spline->y_basis, y, &b, &db, &d2b);

  complex_float* restrict coefs = spline->coefs;
  int xs = spline->x_stride;
  int xs2 = 2*xs;
  #define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  float* restrict p = (float*)P(0,0);
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  p += xs2;
  _mm_prefetch ((const char*)p, _MM_HINT_T0);  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);

  p = (float *)P(0,0);
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (p);    tmp1 = _mm_loadu_ps (p+4); p+= xs2;
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));

  _MM_MATVEC4_PS (r0, r1, r2, r3,   b,   bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   b,   bPi);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  db,  dbPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  db,  dbPi);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2b, d2bPr);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2b, d2bPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;
  float *hess_d2x_r  = ((float*)hess)+0;
  float *hess_d2x_i  = ((float*)hess)+1;
  float *hess_d2y_r  = ((float*)hess)+6;
  float *hess_d2y_i  = ((float*)hess)+7;
  float *hess_dxdy_r = ((float*)hess)+2;
  float *hess_dxdy_i = ((float*)hess)+3;

  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bPr, *gradr0);
  _MM_DOT4_PS (da, bPi, *gradr0);
  _MM_DOT4_PS (a, dbPr, *gradi1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
  // Compute Hessian
  _MM_DOT4_PS (d2a, bPr, *hess_d2x_r);
  _MM_DOT4_PS (d2a, bPi, *hess_d2x_i);
  _MM_DOT4_PS (a, d2bPr, *hess_d2y_r);
  _MM_DOT4_PS (a, d2bPi, *hess_d2y_i);
  _MM_DOT4_PS (da, dbPr, *hess_dxdy_r);
  _MM_DOT4_PS (da, dbPi, *hess_dxdy_i);
#undef P
}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_NUBspline_3d_c (NUBspline_3d_c * restrict spline, 
		    double x, double y, double z,
		    complex_float* restrict val)
{
__m128 a, b, c, cPr[4], bcPr, cPi[4], bcPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;
  int ix = get_NUBasis_funcs_sse_s (spline->x_basis, x, &a);
  int iy = get_NUBasis_funcs_sse_s (spline->y_basis, y, &b);
  int iz = get_NUBasis_funcs_sse_s (spline->z_basis, z, &c);

  int xs  = spline->x_stride;
  int ys  = spline->y_stride;
  int ys2 = 2*ys;
  int ys3 = 3*ys;
  complex_float* restrict coefs = spline->coefs;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  complex_float *p = (complex_float*)P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)p      , _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ ys), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  p = (complex_float*)P(0,0,0);
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  p += xs;
  // 2nd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  p += xs;
  // 3rd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  p += xs;
  // 4th quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],   b,   bcPr);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],   b,   bcPi);


  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  
  // Compute value
  _MM_DOT4_PS (a, bcPr, *valr);
  _MM_DOT4_PS (a, bcPi, *vali);
#undef P
}

/* Value and gradient */
inline void
eval_NUBspline_3d_c_vg (NUBspline_3d_c * restrict spline, 
			double x, double y, double z,
			complex_float* restrict val, complex_float* restrict grad)
{
  __m128 a, b, c, da, db, dc,
    cPr[4], dcPr[4],  bcPr, dbcPr, bdcPr, 
    cPi[4], dcPi[4],  bcPi, dbcPi, bdcPi, 
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  int ix = get_NUBasis_dfuncs_sse_s (spline->x_basis, x, &a, &da);
  int iy = get_NUBasis_dfuncs_sse_s (spline->y_basis, y, &b, &db);
  int iz = get_NUBasis_dfuncs_sse_s (spline->z_basis, z, &c, &dc);

  int xs  = spline->x_stride;
  int ys  = spline->y_stride;
  int ys2 = 2*ys;
  int ys3 = 3*ys;
  complex_float* restrict coefs = spline->coefs;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  complex_float *p = (complex_float*)P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)p      , _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ ys), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  p = (complex_float*)P(0,0,0);
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  p += xs;
  // 2nd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  p += xs;
  // 3rd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  p += xs;
  // 4th quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],   b,   bcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],  db,  dbcPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],   b,  bdcPr);

  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],   b,   bcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],  db,  dbcPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],   b,  bdcPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;
  float *gradr2 = ((float *)grad)+4;
  float *gradi2 = ((float *)grad)+5;
  
  // Compute value
  _MM_DOT4_PS (a, bcPr, *valr);
  _MM_DOT4_PS (a, bcPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bcPr, *gradr0);
  _MM_DOT4_PS (a, dbcPr, *gradr1);
  _MM_DOT4_PS (a, bdcPr, *gradr2);
  _MM_DOT4_PS (da, bcPi, *gradi0);
  _MM_DOT4_PS (a, dbcPi, *gradi1);
  _MM_DOT4_PS (a, bdcPi, *gradi2);
#undef P
}



/* Value, gradient, and laplacian */
inline void
eval_NUBspline_3d_c_vgl (NUBspline_3d_c * restrict spline, 
			 double x, double y, double z,
			 complex_float* restrict val, complex_float* restrict grad, 
			 complex_float* restrict lapl)
{
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cPr[4], dcPr[4], d2cPr[4], bcPr, dbcPr, bdcPr, d2bcPr, bd2cPr,
    cPi[4], dcPi[4], d2cPi[4], bcPi, dbcPi, bdcPi, d2bcPi, bd2cPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  int ix = get_NUBasis_d2funcs_sse_s (spline->x_basis, x, &a, &da, &d2a);
  int iy = get_NUBasis_d2funcs_sse_s (spline->y_basis, y, &b, &db, &d2b);
  int iz = get_NUBasis_d2funcs_sse_s (spline->z_basis, z, &c, &dc, &d2c);

  int xs  = spline->x_stride;
  int ys  = spline->y_stride;
  int ys2 = 2*ys;
  int ys3 = 3*ys;
  complex_float* restrict coefs = spline->coefs;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  complex_float *p = (complex_float*)P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)p      , _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ ys), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  p = (complex_float*)P(0,0,0);
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[0]);
  p += xs;
  // 2nd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[1]);
  p += xs;
  // 3rd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[2]);
  p += xs;
  // 4th quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],   b,   bcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],  db,  dbcPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],   b,  bdcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3], d2b, d2bcPr);
  _MM_MATVEC4_PS (d2cPr[0], d2cPr[1], d2cPr[2], d2cPr[3],   b, bd2cPr);

  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],   b,   bcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],  db,  dbcPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],   b,  bdcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3], d2b, d2bcPi);
  _MM_MATVEC4_PS (d2cPi[0], d2cPi[1], d2cPi[2], d2cPi[3],   b, bd2cPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;
  float *gradr2 = ((float *)grad)+4;
  float *gradi2 = ((float *)grad)+5;
  float d2x_r, d2x_i, d2y_r, d2y_i, d2z_r, d2z_i;
  
  // Compute value
  _MM_DOT4_PS (a, bcPr, *valr);
  _MM_DOT4_PS (a, bcPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bcPr, *gradr0);
  _MM_DOT4_PS (a, dbcPr, *gradr1);
  _MM_DOT4_PS (a, bdcPr, *gradr2);
  _MM_DOT4_PS (da, bcPi, *gradi0);
  _MM_DOT4_PS (a, dbcPi, *gradi1);
  _MM_DOT4_PS (a, bdcPi, *gradi2);
  // Compute laplacian
  _MM_DOT4_PS (d2a, bcPr, d2x_r);
  _MM_DOT4_PS (a, d2bcPr, d2y_r);
  _MM_DOT4_PS (a, bd2cPr, d2z_r);
  _MM_DOT4_PS (d2a, bcPi, d2x_i);
  _MM_DOT4_PS (a, d2bcPi, d2y_i);
  _MM_DOT4_PS (a, bd2cPi, d2z_i);
#ifdef __cplusplus
  *lapl = std::complex<float>(d2x_r + d2y_r + d2z_r,
			      d2x_i + d2y_i + d2z_i);
#else
  *lapl = (d2x_r + d2y_r + d2z_r) + 1.0if * (d2x_i+d2y_i+d2z_i);
#endif
#undef P
}





/* Value, gradient, and Hessian */
inline void
eval_NUBspline_3d_c_vgh (NUBspline_3d_c * restrict spline, 
			 double x, double y, double z,
			 complex_float* restrict val, complex_float* restrict grad, 
			 complex_float* restrict hess)
{
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cPr[4], dcPr[4], d2cPr[4], bcPr, dbcPr, bdcPr, d2bcPr, dbdcPr, bd2cPr,
    cPi[4], dcPi[4], d2cPi[4], bcPi, dbcPi, bdcPi, d2bcPi, dbdcPi, bd2cPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  int ix = get_NUBasis_d2funcs_sse_s (spline->x_basis, x, &a, &da, &d2a);
  int iy = get_NUBasis_d2funcs_sse_s (spline->y_basis, y, &b, &db, &d2b);
  int iz = get_NUBasis_d2funcs_sse_s (spline->z_basis, z, &c, &dc, &d2c);

  int xs  = spline->x_stride;
  int ys  = spline->y_stride;
  int ys2 = 2*ys;
  int ys3 = 3*ys;
  complex_float* restrict coefs = spline->coefs;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  complex_float *p = (complex_float*)P(0,0,0);
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)p      , _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ ys), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);
  p += xs;
  _mm_prefetch ((const char*)(p    ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p    +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys ), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys +2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys2), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys2+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+ys3), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+ys3+2), _MM_HINT_T0);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  p = (complex_float*)P(0,0,0);
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[0]);
  p += xs;
  // 2nd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[1]);
  p += xs;
  // 3rd quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[2]);
  p += xs;
  // 4th quarter
  tmp0 = _mm_loadu_ps ((float*)(p  ));    tmp1 = _mm_loadu_ps ((float*)(p+2));  
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys ));    tmp1 = _mm_loadu_ps ((float*)(p+ys+2));  
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys2 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys2+2));  
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps ((float*)(p+ys3 ));    tmp1 = _mm_loadu_ps ((float*)(p+ys3+2));  
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[3]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[3]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],   b,   bcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3],  db,  dbcPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],   b,  bdcPr);
  _MM_MATVEC4_PS (  cPr[0],   cPr[1],   cPr[2],   cPr[3], d2b, d2bcPr);
  _MM_MATVEC4_PS (d2cPr[0], d2cPr[1], d2cPr[2], d2cPr[3],   b, bd2cPr);
  _MM_MATVEC4_PS ( dcPr[0],  dcPr[1],  dcPr[2],  dcPr[3],  db, dbdcPr);

  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],   b,   bcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3],  db,  dbcPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],   b,  bdcPi);
  _MM_MATVEC4_PS (  cPi[0],   cPi[1],   cPi[2],   cPi[3], d2b, d2bcPi);
  _MM_MATVEC4_PS (d2cPi[0], d2cPi[1], d2cPi[2], d2cPi[3],   b, bd2cPi);
  _MM_MATVEC4_PS ( dcPi[0],  dcPi[1],  dcPi[2],  dcPi[3],  db, dbdcPi);

  float *valr   = ((float*)val)  +0;
  float *vali   = ((float*)val)  +1;
  float *gradr0 = ((float *)grad)+0;
  float *gradi0 = ((float *)grad)+1;
  float *gradr1 = ((float *)grad)+2;
  float *gradi1 = ((float *)grad)+3;
  float *gradr2 = ((float *)grad)+4;
  float *gradi2 = ((float *)grad)+5;
  
  // Compute value
  _MM_DOT4_PS (a, bcPr, *valr);
  _MM_DOT4_PS (a, bcPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bcPr, *gradr0);
  _MM_DOT4_PS (a, dbcPr, *gradr1);
  _MM_DOT4_PS (a, bdcPr, *gradr2);
  _MM_DOT4_PS (da, bcPi, *gradi0);
  _MM_DOT4_PS (a, dbcPi, *gradi1);
  _MM_DOT4_PS (a, bdcPi, *gradi2);
  // Compute hessian
  _MM_DOT4_PS (d2a, bcPr, *(float*)(&hess[0]));
  _MM_DOT4_PS (a, d2bcPr, *(float*)(&hess[4]));
  _MM_DOT4_PS (a, bd2cPr, *(float*)(&hess[8]));
  _MM_DOT4_PS (da, dbcPr, *(float*)(&hess[1]));
  _MM_DOT4_PS (da, bdcPr, *(float*)(&hess[2]));
  _MM_DOT4_PS (a, dbdcPr, *(float*)(&hess[5]));

  _MM_DOT4_PS (d2a, bcPi, *((float*)(&hess[0])+1));
  _MM_DOT4_PS (a, d2bcPi, *((float*)(&hess[4])+1));
  _MM_DOT4_PS (a, bd2cPi, *((float*)(&hess[8])+1));
  _MM_DOT4_PS (da, dbcPi, *((float*)(&hess[1])+1));
  _MM_DOT4_PS (da, bdcPi, *((float*)(&hess[2])+1));
  _MM_DOT4_PS (a, dbdcPi, *((float*)(&hess[5])+1));


  // Multiply gradients and hessians by appropriate grid inverses
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];

#undef P
}
#undef _MM_MATVEC4_PS
#undef _MM_DOT4_PS
#endif
