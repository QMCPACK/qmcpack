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

#ifndef NUBSPLINE_EVAL_SSE_Z_H
#define NUBSPLINE_EVAL_SSE_Z_H

#include <math.h>
#include <stdio.h>
#include "nubspline_structs.h"

#ifdef HAVE_SSE2
#include <xmmintrin.h>
#include <emmintrin.h>
#endif


#ifdef HAVE_SSE3
#include <pmmintrin.h>
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
inline void
eval_NUBspline_1d_z (NUBspline_1d_z * restrict spline, 
		     double x, complex_double* restrict val)
{
  double bfuncs[4];
  int i = get_NUBasis_funcs_d (spline->x_basis, x, bfuncs);
  complex_double* restrict coefs = spline->coefs;
  *val = (coefs[i+0]*bfuncs[0] +coefs[i+1]*bfuncs[1] +
	  coefs[i+2]*bfuncs[2] +coefs[i+3]*bfuncs[3]);
}

/* Value and first derivative */
inline void
eval_NUBspline_1d_z_vg (NUBspline_1d_z * restrict spline, double x, 
			complex_double* restrict val, complex_double* restrict grad)
{
  double bfuncs[4], dbfuncs[4];
  int i = get_NUBasis_dfuncs_d (spline->x_basis, x, bfuncs, dbfuncs);
  complex_double* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]* bfuncs[0] + coefs[i+1]* bfuncs[1] +
	   coefs[i+2]* bfuncs[2] + coefs[i+3]* bfuncs[3]);
  *grad = (coefs[i+0]*dbfuncs[0] + coefs[i+1]*dbfuncs[1] +
	   coefs[i+2]*dbfuncs[2] + coefs[i+3]*dbfuncs[3]);
}

/* Value, first derivative, and second derivative */
inline void
eval_NUBspline_1d_z_vgl (NUBspline_1d_z * restrict spline, double x, 
			complex_double* restrict val, complex_double* restrict grad,
			complex_double* restrict lapl)
{
  double bfuncs[4], dbfuncs[4], d2bfuncs[4];
  int i = get_NUBasis_d2funcs_d (spline->x_basis, x, bfuncs, dbfuncs, d2bfuncs);
  complex_double* restrict coefs = spline->coefs;
  *val =  (coefs[i+0]*  bfuncs[0] + coefs[i+1]*  bfuncs[1] +
	   coefs[i+2]*  bfuncs[2] + coefs[i+3]*  bfuncs[3]);
  *grad = (coefs[i+0]* dbfuncs[0] + coefs[i+1]* dbfuncs[1] +
	   coefs[i+2]* dbfuncs[2] + coefs[i+3]* dbfuncs[3]);
  *lapl = (coefs[i+0]*d2bfuncs[0] + coefs[i+1]*d2bfuncs[1] +
	   coefs[i+2]*d2bfuncs[2] + coefs[i+3]*d2bfuncs[3]);

}

inline void
eval_NUBspline_1d_z_vgh (NUBspline_1d_z * restrict spline, double x, 
			complex_double* restrict val, complex_double* restrict grad,
			complex_double* restrict hess)
{
  eval_NUBspline_1d_z_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_NUBspline_2d_z (NUBspline_2d_z * restrict spline, 
		    double x, double y, complex_double* restrict val)
{
  __m128d
    a01, b01, a23, b23, bP01r, bP23r, bP01i, bP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
  
  int ix = get_NUBasis_funcs_sse_d (spline->x_basis, x, &a01, &a23);
  int iy = get_NUBasis_funcs_sse_d (spline->y_basis, y, &b01, &b23);
  int xs   = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*)P(0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;

  tmp0 = _mm_load_pd (P(0,0));  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);

  tmp0 = _mm_load_pd (P(2,0));  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));  tmp1 = _mm_load_pd (P(3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bP01r, bP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bP01i, bP23i, *((double*)val+1));

#undef P
}


/* Value and gradient */
inline void
eval_NUBspline_2d_z_vg (NUBspline_2d_z * restrict spline, 
		       double x, double y, 
		       complex_double* restrict val, complex_double* restrict grad)
{
  __m128d
    a01, b01, da01, db01, 
    a23, b23, da23, db23, 
    bP01r, dbP01r, bP23r, dbP23r, 
    bP01i, dbP01i, bP23i, dbP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
  
  int ix = get_NUBasis_dfuncs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23);
  int iy = get_NUBasis_dfuncs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23);
  int xs   = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*)P(0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;

  tmp0 = _mm_load_pd (P(0,0));  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);

  tmp0 = _mm_load_pd (P(2,0));  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));  tmp1 = _mm_load_pd (P(3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP23i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bP01r, bP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bP01i, bP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bP01r,  bP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bP01i,  bP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbP01r, dbP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbP01i, dbP23i, dgrad[3]);

#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_NUBspline_2d_z_vgl (NUBspline_2d_z * restrict spline, 
			double x, double y, complex_double* restrict val, 
			complex_double* restrict grad, complex_double* restrict lapl)
{
 __m128d
    a01, b01, da01, db01, d2a01, d2b01,
    a23, b23, da23, db23, d2a23, d2b23,
    bP01r, dbP01r, d2bP01r, 
    bP23r, dbP23r, d2bP23r, 
    bP01i, dbP01i, d2bP01i, 
    bP23i, dbP23i, d2bP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
  
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int xs   = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*)P(0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;

  tmp0 = _mm_load_pd (P(0,0));  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01i);

  tmp0 = _mm_load_pd (P(2,0));  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));  tmp1 = _mm_load_pd (P(3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP23i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP23i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bP01r, bP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bP01i, bP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bP01r,  bP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bP01i,  bP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbP01r, dbP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbP01i, dbP23i, dgrad[3]);
  // Compute Laplacian
  double d2x_r, d2x_i, d2y_r, d2y_i;
  _MM_DOT4_PD (d2a01, d2a23, bP01r, bP23r, d2x_r);
  _MM_DOT4_PD (d2a01, d2a23, bP01i, bP23i, d2x_i);
  _MM_DOT4_PD (a01, a23, d2bP01r, d2bP23r, d2y_r);
  _MM_DOT4_PD (a01, a23, d2bP01i, d2bP23i, d2y_i);
#ifdef __cplusplus
  *lapl = std::complex<double>(d2x_r + d2y_r,
			       d2x_i + d2y_i);
#else
  *lapl = (d2x_r + d2y_r) + 1.0i*(d2x_i + d2y_i);
#endif
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_NUBspline_2d_z_vgh (NUBspline_2d_z * restrict spline, 
			double x, double y, complex_double* restrict val, 
			complex_double* restrict grad, complex_double* restrict hess)
{
  __m128d
    a01, b01, da01, db01, d2a01, d2b01,
    a23, b23, da23, db23, d2a23, d2b23,
    bP01r, dbP01r, d2bP01r, 
    bP23r, dbP23r, d2bP23r, 
    bP01i, dbP01i, d2bP01i, 
    bP23i, dbP23i, d2bP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
  
  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int xs   = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*)P(0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);  _mm_prefetch ((const char*)(p+6), _MM_HINT_T0); p+= xs;

  tmp0 = _mm_load_pd (P(0,0));  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01i);

  tmp0 = _mm_load_pd (P(2,0));  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));  tmp1 = _mm_load_pd (P(3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP23i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP23i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP23r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bP01r, bP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bP01i, bP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  double *dhess = (double*) hess;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bP01r,  bP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bP01i,  bP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbP01r, dbP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbP01i, dbP23i, dgrad[3]);
  // Compute Hessian
  _MM_DOT4_PD (d2a01, d2a23, bP01r, bP23r, dhess[0]);
  _MM_DOT4_PD (d2a01, d2a23, bP01i, bP23i, dhess[1]);
  _MM_DOT4_PD (a01, a23, d2bP01r, d2bP23r, dhess[6]);
  _MM_DOT4_PD (a01, a23, d2bP01i, d2bP23i, dhess[7]);
  _MM_DOT4_PD (da01, da23, dbP01r, dbP23r, dhess[2]);
  _MM_DOT4_PD (da01, da23, dbP01i, dbP23i, dhess[3]);
  
  hess[2] = hess[1];
#undef P
}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_NUBspline_3d_z (NUBspline_3d_z * restrict spline, 
		    double x, double y, double z,
		    complex_double* restrict val)
{
 __m128d 
    a01, b01, c01, a23, b23, c23, cPr[8], cPi[8], bcP01r, bcP23r, bcP01i, bcP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;

  int ix = get_NUBasis_funcs_sse_d (spline->x_basis, x, &a01, &a23);
  int iy = get_NUBasis_funcs_sse_d (spline->y_basis, y, &b01, &b23);
  int iz = get_NUBasis_funcs_sse_d (spline->z_basis, z, &c01, &c23);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
  int xs2 = 2*xs;
  int ys2 = 2*ys;
  int delta = xs2-3*ys2;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*) P(0,0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  

  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  p = (double*) P(0,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  
  // 2nd eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  
  // 3rd eighth
  p = (double*) P(1,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);

  // 4th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);

  // 5th eighth
  p = (double*) P(2,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);

  // 6th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);

  // 7th eighth
  p = (double*) P(3,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);

  // 8th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6); p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[7]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[0], cPr[1], cPr[2], cPr[3], bcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[0], cPi[1], cPi[2], cPi[3], bcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[4], cPr[5], cPr[6], cPr[7], bcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[4], cPi[5], cPi[6], cPi[7], bcP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01r, bcP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bcP01i, bcP23i, *((double*)val+1));

#undef P
}

/* Value and gradient */
inline void
eval_NUBspline_3d_z_vg (NUBspline_3d_z * restrict spline, 
			double x, double y, double z,
			complex_double* restrict val, complex_double* restrict grad)
{
  __m128d 
    a01, b01, c01, da01, db01, dc01, 
    a23, b23, c23, da23, db23, dc23, 
    cPr[8], dcPr[8], 
    cPi[8], dcPi[8], 
    bcP01r, dbcP01r, bdcP01r, 
    bcP23r, dbcP23r, bdcP23r, 
    bcP01i, dbcP01i, bdcP01i, 
    bcP23i, dbcP23i, bdcP23i, 
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;

  int ix = get_NUBasis_dfuncs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23);
  int iy = get_NUBasis_dfuncs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23);
  int iz = get_NUBasis_dfuncs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
  int xs2 = 2*xs;
  int ys2 = 2*ys;
  int delta = xs2-3*ys2;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*) P(0,0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  

  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  p = (double*) P(0,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  
  // 2nd eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  
  // 3rd eighth
  p = (double*) P(1,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);

  // 4th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);

  // 5th eighth
  p = (double*) P(2,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);

  // 6th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);

  // 7th eighth
  p = (double*) P(3,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);

  // 8th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6); p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[7]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[7]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[0], cPr[1], cPr[2], cPr[3], bcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[0], cPi[1], cPi[2], cPi[3], bcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[4], cPr[5], cPr[6], cPr[7], bcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[4], cPi[5], cPi[6], cPi[7], bcP23i);

  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[0], cPr[1], cPr[2], cPr[3], dbcP01r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[0], cPi[1], cPi[2], cPi[3], dbcP01i);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[4], cPr[5], cPr[6], cPr[7], dbcP23r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[4], cPi[5], cPi[6], cPi[7], dbcP23i);

  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[0], dcPr[1], dcPr[2], dcPr[3], bdcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[0], dcPi[1], dcPi[2], dcPi[3], bdcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[4], dcPr[5], dcPr[6], dcPr[7], bdcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[4], dcPi[5], dcPi[6], dcPi[7], bdcP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01r, bcP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bcP01i, bcP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bcP01r,  bcP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bcP01i,  bcP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbcP01r, dbcP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbcP01i, dbcP23i, dgrad[3]);
  _MM_DOT4_PD ( a01,  a23, bdcP01r, bdcP23r, dgrad[4]);
  _MM_DOT4_PD ( a01,  a23, bdcP01i, bdcP23i, dgrad[5]);

#undef P

}



/* Value, gradient, and laplacian */
inline void
eval_NUBspline_3d_z_vgl (NUBspline_3d_z * restrict spline, 
			 double x, double y, double z,
			 complex_double* restrict val, complex_double* restrict grad, 
			 complex_double* restrict lapl)
{
  __m128d 
    a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
    a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
    cPr[8], dcPr[8], d2cPr[8], 
    cPi[8], dcPi[8], d2cPi[8], 
    bcP01r, dbcP01r, bdcP01r, d2bcP01r, bd2cP01r,
    bcP23r, dbcP23r, bdcP23r, d2bcP23r, bd2cP23r,
    bcP01i, dbcP01i, bdcP01i, d2bcP01i, bd2cP01i,
    bcP23i, dbcP23i, bdcP23i, d2bcP23i, bd2cP23i,
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;

  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int iz = get_NUBasis_d2funcs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23, &d2c01, &d2c23);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
  int xs2 = 2*xs;
  int ys2 = 2*ys;
  int delta = xs2-3*ys2;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*) P(0,0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  

  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  p = (double*) P(0,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[0]);
  
  // 2nd eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[1]);
  
  // 3rd eighth
  p = (double*) P(1,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[2]);

  // 4th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[3]);

  // 5th eighth
  p = (double*) P(2,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[4]);

  // 6th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[5]);

  // 7th eighth
  p = (double*) P(3,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[6]);

  // 8th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6); p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[7]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[7]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[7]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[0], cPr[1], cPr[2], cPr[3], bcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[0], cPi[1], cPi[2], cPi[3], bcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[4], cPr[5], cPr[6], cPr[7], bcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[4], cPi[5], cPi[6], cPi[7], bcP23i);

  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[0], cPr[1], cPr[2], cPr[3], dbcP01r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[0], cPi[1], cPi[2], cPi[3], dbcP01i);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[4], cPr[5], cPr[6], cPr[7], dbcP23r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[4], cPi[5], cPi[6], cPi[7], dbcP23i);

  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[0], dcPr[1], dcPr[2], dcPr[3], bdcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[0], dcPi[1], dcPi[2], dcPi[3], bdcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[4], dcPr[5], dcPr[6], dcPr[7], bdcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[4], dcPi[5], dcPi[6], dcPi[7], bdcP23i);

  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPr[0], cPr[1], cPr[2], cPr[3], d2bcP01r);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPi[0], cPi[1], cPi[2], cPi[3], d2bcP01i);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPr[4], cPr[5], cPr[6], cPr[7], d2bcP23r);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPi[4], cPi[5], cPi[6], cPi[7], d2bcP23i);

  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPr[0], d2cPr[1], d2cPr[2], d2cPr[3], bd2cP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPi[0], d2cPi[1], d2cPi[2], d2cPi[3], bd2cP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPr[4], d2cPr[5], d2cPr[6], d2cPr[7], bd2cP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPi[4], d2cPi[5], d2cPi[6], d2cPi[7], bd2cP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01r, bcP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bcP01i, bcP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bcP01r,  bcP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bcP01i,  bcP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbcP01r, dbcP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbcP01i, dbcP23i, dgrad[3]);
  _MM_DOT4_PD ( a01,  a23, bdcP01r, bdcP23r, dgrad[4]);
  _MM_DOT4_PD ( a01,  a23, bdcP01i, bdcP23i, dgrad[5]);
  
  // Compute Laplacian
  double d2x_r, d2x_i, d2y_r, d2y_i, d2z_r, d2z_i;
  // d2x
  _MM_DOT4_PD (d2a01, d2a23, bcP01r, bcP23r, d2x_r);
  _MM_DOT4_PD (d2a01, d2a23, bcP01i, bcP23i, d2x_i);
  // d2y
  _MM_DOT4_PD (a01, a23, d2bcP01r, d2bcP23r, d2y_r);
  _MM_DOT4_PD (a01, a23, d2bcP01i, d2bcP23i, d2y_i);
  // d2z
  _MM_DOT4_PD (a01, a23, bd2cP01r, bd2cP23r, d2z_r);
  _MM_DOT4_PD (a01, a23, bd2cP01i, bd2cP23i, d2z_i);
#ifdef __cplusplus
  *lapl = std::complex<double>(d2x_r + d2y_r + d2z_r,
			       d2x_i + d2y_i + d2z_i);
#else
  *lapl = (d2x_r + d2y_r + d2z_r) + 1.0i*(d2x_i + d2y_i + d2z_i);
#endif
#undef P

}





/* Value, gradient, and Hessian */
inline void
eval_NUBspline_3d_z_vgh (NUBspline_3d_z * restrict spline, 
			 double x, double y, double z,
			 complex_double* restrict val, complex_double* restrict grad, 
			 complex_double* restrict hess)
{
  __m128d 
    a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
    a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
   cPr[8], dcPr[8], d2cPr[8], 
    cPi[8], dcPi[8], d2cPi[8], 
    bcP01r, dbcP01r, bdcP01r, d2bcP01r, dbdcP01r, bd2cP01r,
    bcP23r, dbcP23r, bdcP23r, d2bcP23r, dbdcP23r, bd2cP23r,
    bcP01i, dbcP01i, bdcP01i, d2bcP01i, dbdcP01i, bd2cP01i,
    bcP23i, dbcP23i, bdcP23i, d2bcP23i, dbdcP23i, bd2cP23i,
    tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;

  int ix = get_NUBasis_d2funcs_sse_d (spline->x_basis, x, &a01, &a23, &da01, &da23, &d2a01, &d2a23);
  int iy = get_NUBasis_d2funcs_sse_d (spline->y_basis, y, &b01, &b23, &db01, &db23, &d2b01, &d2b23);
  int iz = get_NUBasis_d2funcs_sse_d (spline->z_basis, z, &c01, &c23, &dc01, &dc23, &d2c01, &d2c23);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
  int xs2 = 2*xs;
  int ys2 = 2*ys;
  int delta = xs2-3*ys2;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  double *restrict p = (double*) P(0,0,0);
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += delta;

  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  p += ys2;
  _mm_prefetch ((const char*)(p+0), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+2), _MM_HINT_T0);
  _mm_prefetch ((const char*)(p+4), _MM_HINT_T0);    _mm_prefetch ((const char*)(p+6), _MM_HINT_T0);  

  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  p = (double*) P(0,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[0]);
  
  // 2nd eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[1]);
  
  // 3rd eighth
  p = (double*) P(1,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[2]);

  // 4th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[3]);

  // 5th eighth
  p = (double*) P(2,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[4]);

  // 6th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[5]);

  // 7th eighth
  p = (double*) P(3,0,0);
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += ys2;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[6]);

  // 8th eighth
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6); p += ys2;
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+0);  tmp1 = _mm_load_pd (p+2);
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (p+4);  tmp1 = _mm_load_pd (p+6);  p += delta;
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[7]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[7]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[7]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[7]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[0], cPr[1], cPr[2], cPr[3], bcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[0], cPi[1], cPi[2], cPi[3], bcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPr[4], cPr[5], cPr[6], cPr[7], bcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, cPi[4], cPi[5], cPi[6], cPi[7], bcP23i);

  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[0], cPr[1], cPr[2], cPr[3], dbcP01r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[0], cPi[1], cPi[2], cPi[3], dbcP01i);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPr[4], cPr[5], cPr[6], cPr[7], dbcP23r);
  _MM_DDOT4_PD (db01, db23, db01, db23, cPi[4], cPi[5], cPi[6], cPi[7], dbcP23i);

  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[0], dcPr[1], dcPr[2], dcPr[3], bdcP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[0], dcPi[1], dcPi[2], dcPi[3], bdcP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPr[4], dcPr[5], dcPr[6], dcPr[7], bdcP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, dcPi[4], dcPi[5], dcPi[6], dcPi[7], bdcP23i);

  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPr[0], cPr[1], cPr[2], cPr[3], d2bcP01r);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPi[0], cPi[1], cPi[2], cPi[3], d2bcP01i);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPr[4], cPr[5], cPr[6], cPr[7], d2bcP23r);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23, cPi[4], cPi[5], cPi[6], cPi[7], d2bcP23i);

  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPr[0], d2cPr[1], d2cPr[2], d2cPr[3], bd2cP01r);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPi[0], d2cPi[1], d2cPi[2], d2cPi[3], bd2cP01i);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPr[4], d2cPr[5], d2cPr[6], d2cPr[7], bd2cP23r);
  _MM_DDOT4_PD (b01, b23, b01, b23, d2cPi[4], d2cPi[5], d2cPi[6], d2cPi[7], bd2cP23i);
  
  _MM_DDOT4_PD (db01, db23, db01, db23, dcPr[0], dcPr[1], dcPr[2], dcPr[3], dbdcP01r);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcPi[0], dcPi[1], dcPi[2], dcPi[3], dbdcP01i);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcPr[4], dcPr[5], dcPr[6], dcPr[7], dbdcP23r);
  _MM_DDOT4_PD (db01, db23, db01, db23, dcPi[4], dcPi[5], dcPi[6], dcPi[7], dbdcP23i);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01r, bcP23r, *((double*)val+0));
  _MM_DOT4_PD (a01, a23, bcP01i, bcP23i, *((double*)val+1));

  double *dgrad = (double*) grad;
  // Compute gradient
  _MM_DOT4_PD (da01, da23,  bcP01r,  bcP23r, dgrad[0]);
  _MM_DOT4_PD (da01, da23,  bcP01i,  bcP23i, dgrad[1]);
  _MM_DOT4_PD ( a01,  a23, dbcP01r, dbcP23r, dgrad[2]);
  _MM_DOT4_PD ( a01,  a23, dbcP01i, dbcP23i, dgrad[3]);
  _MM_DOT4_PD ( a01,  a23, bdcP01r, bdcP23r, dgrad[4]);
  _MM_DOT4_PD ( a01,  a23, bdcP01i, bdcP23i, dgrad[5]);
  
  double *dhess = (double*) hess;
  // Compute hessian
  // d2x
  _MM_DOT4_PD (d2a01, d2a23, bcP01r, bcP23r, dhess[0]);
  _MM_DOT4_PD (d2a01, d2a23, bcP01i, bcP23i, dhess[1]);
  // d2y
  _MM_DOT4_PD (a01, a23, d2bcP01r, d2bcP23r, dhess[8]);
  _MM_DOT4_PD (a01, a23, d2bcP01i, d2bcP23i, dhess[9]);
  // d2z
  _MM_DOT4_PD (a01, a23, bd2cP01r, bd2cP23r, dhess[16]);
  _MM_DOT4_PD (a01, a23, bd2cP01i, bd2cP23i, dhess[17]);
  // dx dy
  _MM_DOT4_PD (da01, da23, dbcP01r, dbcP23r, dhess[2]);
  _MM_DOT4_PD (da01, da23, dbcP01i, dbcP23i, dhess[3]);
  // dx dz
  _MM_DOT4_PD (da01, da23, bdcP01r, bdcP23r, dhess[4]);
  _MM_DOT4_PD (da01, da23, bdcP01i, bdcP23i, dhess[5]);
  // dy dz
  _MM_DOT4_PD (a01, a23, dbdcP01r, dbdcP23r, dhess[10]);
  _MM_DOT4_PD (a01, a23, dbdcP01i, dbdcP23i, dhess[11]);
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];

#undef P

}

#undef _MM_DDOT4_PD
#undef _MM_DOT4_PD

#endif
