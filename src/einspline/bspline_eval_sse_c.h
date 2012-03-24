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

#ifndef BSPLINE_EVAL_SSE_C_H
#define BSPLINE_EVAL_SSE_C_H

#include "bspline_structs.h"

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
  #include <pmmintrin.h>
#endif
#include <math.h>

// extern __m128   A0,   A1,   A2,   A3;
// extern __m128  dA0,  dA1,  dA2,  dA3;
// extern __m128 d2A0, d2A1, d2A2, d2A3;
extern __m128* restrict A_s;


inline void
print__m128 (__m128 val)
{
  float v[4];
  __m128 vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(0,0,0,0));
  _mm_store_ss (&(v[0]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(1,1,1,1));
  _mm_store_ss (&(v[1]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(2,2,2,2));
  _mm_store_ss (&(v[2]), vshuf);
  vshuf = _mm_shuffle_ps (val, val, _MM_SHUFFLE(3,3,3,3));
  _mm_store_ss (&(v[3]), vshuf);
  
  fprintf (stderr, "[ %8.5f, %8.5f, %8.5f, %8.5f ]", v[0], v[1], v[2], v[3]);
}



/// SSE3 adds "horizontal add" instructions, which makes things
/// simpler and faster
#ifdef HAVE_SSE3
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                         \
do {                                                                 \
  __m128 _r0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
  __m128 _r1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
  r = _mm_hadd_ps (_r0, _r1);                                        \
 } while (0);
#define _MM_DOT4_PS(_A, _B, _p)                                      \
do {                                                                 \
  __m128 t  = _mm_mul_ps (_A, _B);                                   \
  __m128 t1 = _mm_hadd_ps (t,t);                                     \
  __m128 r  = _mm_hadd_ps (t1, t1);                                  \
  _mm_store_ss (&(_p), r);                                           \
} while(0);
#else
// Use plain-old SSE instructions
#define _MM_MATVEC4_PS(_M0, _M1, _M2, _M3, _v, _r)                   \
do {                                                                 \
  __m128 _r0 = _mm_mul_ps (_M0, _v);                                 \
  __m128 _r1 = _mm_mul_ps (_M1, _v);				     \
  __m128 _r2 = _mm_mul_ps (_M2, _v);                                 \
  __m128 _r3 = _mm_mul_ps (_M3, _v);				     \
  _MM_TRANSPOSE4_PS (_r0, _r1, _r2, _r3);                            \
  _r = _mm_add_ps (_mm_add_ps (_r0, _r1), _mm_add_ps (_r2, _r3));    \
 } while (0);
#define _MM_DOT4_PS(_A, _B, _p)                                      \
do {                                                                 \
  __m128 _t   = _mm_mul_ps (_A, _B);                                 \
  __m128 alo  = _mm_shuffle_ps (_t, _t, _MM_SHUFFLE(0,1,0,1));	     \
  __m128 ahi  = _mm_shuffle_ps (_t, _t, _MM_SHUFFLE(2,3,2,3));	     \
  __m128 _a   = _mm_add_ps (alo, ahi);                               \
  __m128 rlo  = _mm_shuffle_ps (_a, _a, _MM_SHUFFLE(0,0,0,0));	     \
  __m128 rhi  = _mm_shuffle_ps (_a, _a, _MM_SHUFFLE(1,1,1,1));	     \
  __m128 r    = _mm_add_ps (rlo, rhi);                               \
  _mm_store_ss (&(_p), r);                                           \
} while(0);
#endif


/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_c (UBspline_1d_c * restrict spline, 
		    double x, complex_float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_c_vg (UBspline_1d_c * restrict spline, double x, 
		     complex_float* restrict val, complex_float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = (float)spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}

/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_c_vgl (UBspline_1d_c* restrict spline, double x, 
			complex_float* restrict val, 
			complex_float* restrict grad,
			complex_float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = (float)spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = (float)(spline->x_grid.delta_inv * spline->x_grid.delta_inv) * 
    (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
     coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
     coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
     coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

inline void
eval_UBspline_1d_c_vgh (UBspline_1d_c* restrict spline, double x, 
			complex_float* restrict val, 
			complex_float* restrict grad,
			complex_float* restrict hess)
{
  eval_UBspline_1d_c_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_c (UBspline_2d_c * restrict spline, 
		    double x, double y, complex_float* restrict val)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  __m128 tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  __m128 tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);

//   /// SSE mesh point determination
//   __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
//   __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
//   __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
//   xy = _mm_sub_ps (xy, x0y0);
//   // ux = (x - x0)/delta_x and same for y
//   __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
//   // intpart = trunc (ux, uy)
//   __m128i intpart  = _mm_cvttps_epi32(uxuy);
//   __m128i ixiy;
//   _mm_storeu_si128 (&ixiy, intpart);
//   // Store to memory for use in C expressions
//   // xmm registers are stored to memory in reverse order
//   int ix = ((int *)&ixiy)[3];
//   int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
//   __m128 ipart  = _mm_cvtepi32_ps (intpart);
//   __m128 txty   = _mm_sub_ps (uxuy, ipart);
//   __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
//   __m128 t2     = _mm_mul_ps (txty, txty);
//   __m128 t3     = _mm_mul_ps (t2, txty);
//   __m128 tpx    = t3;
//   __m128 tpy    = t2;
//   __m128 tpz    = txty;
//   __m128 zero   = one;
//   _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, bPr, bPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);

  tmp0 = _mm_loadu_ps (P(0,0));  
  tmp1 = _mm_loadu_ps (P(0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,0));  
  tmp1 = _mm_loadu_ps (P(1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,0));  
  tmp1 = _mm_loadu_ps (P(2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,0));  
  tmp1 = _mm_loadu_ps (P(3,2));
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
eval_UBspline_2d_c_vg (UBspline_2d_c * restrict spline, 
		       double x, double y, 
		       complex_float* restrict val, complex_float* restrict grad)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  __m128 tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  __m128 tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);

//   /// SSE mesh point determination
//   __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
//   __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
//   __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
//   xy = _mm_sub_ps (xy, x0y0);
//   // ux = (x - x0)/delta_x and same for y
//   __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
//   // intpart = trunc (ux, uy)
//   __m128i intpart  = _mm_cvttps_epi32(uxuy);
//   __m128i ixiy;
//   _mm_storeu_si128 (&ixiy, intpart);
//   // Store to memory for use in C expressions
//   // xmm registers are stored to memory in reverse order
//   int ix = ((int *)&ixiy)[3];
//   int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
//   __m128 ipart  = _mm_cvtepi32_ps (intpart);
//   __m128 txty   = _mm_sub_ps (uxuy, ipart);
//   __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
//   __m128 t2     = _mm_mul_ps (txty, txty);
//   __m128 t3     = _mm_mul_ps (t2, txty);
//   __m128 tpx    = t3;
//   __m128 tpy    = t2;
//   __m128 tpz    = txty;
//   __m128 zero   = one;
//   _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, da, db, bPr, dbPr, bPi, dbPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpx,  da);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpy,  db);

  tmp0 = _mm_loadu_ps (P(0,0));  
  tmp1 = _mm_loadu_ps (P(0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,0));  
  tmp1 = _mm_loadu_ps (P(1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,0));  
  tmp1 = _mm_loadu_ps (P(2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,0));  
  tmp1 = _mm_loadu_ps (P(3,2));
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
  _MM_DOT4_PS (da, bPi, *gradi0);
  _MM_DOT4_PS (a, dbPr, *gradr1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
  
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;

#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_c_vgl (UBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, complex_float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  __m128 tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  __m128 tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);

//   /// SSE mesh point determination
//   __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
//   __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
//   __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
//   xy = _mm_sub_ps (xy, x0y0);
//   // ux = (x - x0)/delta_x and same for y
//   __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
//   // intpart = trunc (ux, uy)
//   __m128i intpart  = _mm_cvttps_epi32(uxuy);
//   __m128i ixiy;
//   _mm_storeu_si128 (&ixiy, intpart);
//   // Store to memory for use in C expressions
//   // xmm registers are stored to memory in reverse order
//   int ix = ((int *)&ixiy)[3];
//   int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
//   __m128 ipart  = _mm_cvtepi32_ps (intpart);
//   __m128 txty   = _mm_sub_ps (uxuy, ipart);
//   __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
//   __m128 t2     = _mm_mul_ps (txty, txty);
//   __m128 t3     = _mm_mul_ps (t2, txty);
//   __m128 tpx    = t3;
//   __m128 tpy    = t2;
//   __m128 tpz    = txty;
//   __m128 zero   = one;
//   _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, da, db, d2a, d2b, 
    bPr, dbPr, d2bPr, bPi, dbPi, d2bPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b);

  tmp0 = _mm_loadu_ps (P(0,0));  
  tmp1 = _mm_loadu_ps (P(0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,0));  
  tmp1 = _mm_loadu_ps (P(1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,0));  
  tmp1 = _mm_loadu_ps (P(2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,0));  
  tmp1 = _mm_loadu_ps (P(3,2));
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
  float  hess_d2x_r, hess_d2x_i, 
    hess_d2y_r, hess_d2y_i;

  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bPr, *gradr0);
  _MM_DOT4_PS (da, bPi, *gradi0);
  _MM_DOT4_PS (a, dbPr, *gradr1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
  // Compute Hessian
  _MM_DOT4_PS (d2a, bPr, hess_d2x_r);
  _MM_DOT4_PS (d2a, bPi, hess_d2x_i);
  _MM_DOT4_PS (a, d2bPr, hess_d2y_r);
  _MM_DOT4_PS (a, d2bPi, hess_d2y_i);
  
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  hess_d2x_r *= dxInv*dxInv;
  hess_d2x_i *= dxInv*dxInv;
  hess_d2y_r *= dyInv*dyInv;
  hess_d2y_i *= dyInv*dyInv;
#ifdef __cplusplus
  *lapl = std::complex<float>(hess_d2x_r + hess_d2y_r, hess_d2x_i + hess_d2y_i);
#else
  *lapl = (hess_d2x_r + hess_d2y_r) + 1.0fI* (hess_d2x_i + hess_d2y_i);
#endif
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_c_vgh (UBspline_2d_c * restrict spline, 
			double x, double y, complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
//   /// SSE mesh point determination
//   __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
//   __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
//   __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
//   xy = _mm_sub_ps (xy, x0y0);
//   // ux = (x - x0)/delta_x and same for y
//   __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
//   // intpart = trunc (ux, uy)
//   __m128i intpart  = _mm_cvttps_epi32(uxuy);
//   __m128i ixiy;
//   _mm_storeu_si128 (&ixiy, intpart);
//   // Store to memory for use in C expressions
//   // xmm registers are stored to memory in reverse order
//   int ix = ((int *)&ixiy)[3];
//   int iy = ((int *)&ixiy)[2];

//   __m128 ipart  = _mm_cvtepi32_ps (intpart);
//   __m128 txty   = _mm_sub_ps (uxuy, ipart);
//   __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
//   __m128 t2     = _mm_mul_ps (txty, txty);
//   __m128 t3     = _mm_mul_ps (t2, txty);
//   __m128 tpx    = t3;
//   __m128 tpy    = t2;
//   __m128 tpz    = txty;
//   __m128 zero   = one;
//   _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i,j) (const float*)(spline->coefs+(ix+(i))*xs+iy+j)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]

  __m128 tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  __m128 tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);
  
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, da, db, d2a, d2b, 
    bPr, dbPr, d2bPr, bPi, dbPi, d2bPi,
    r0, r1, r2, r3, i0, i1, i2, i3, tmp0, tmp1;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b);

  tmp0 = _mm_loadu_ps (P(0,0));  
  tmp1 = _mm_loadu_ps (P(0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,0));  
  tmp1 = _mm_loadu_ps (P(1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,0));  
  tmp1 = _mm_loadu_ps (P(2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,0));  
  tmp1 = _mm_loadu_ps (P(3,2));
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
  float *hess_dydx_r = ((float*)hess)+4;
  float *hess_dydx_i = ((float*)hess)+5;
  // Compute value
  _MM_DOT4_PS (a, bPr, *valr);
  _MM_DOT4_PS (a, bPi, *vali);
  // Compute gradient
  _MM_DOT4_PS (da, bPr, *gradr0);
  _MM_DOT4_PS (da, bPi, *gradi0);
  _MM_DOT4_PS (a, dbPr, *gradr1);
  _MM_DOT4_PS (a, dbPi, *gradi1);
  // Compute Hessian
  _MM_DOT4_PS (d2a, bPr, *hess_d2x_r);
  _MM_DOT4_PS (d2a, bPi, *hess_d2x_i);
  _MM_DOT4_PS (a, d2bPr, *hess_d2y_r);
  _MM_DOT4_PS (a, d2bPi, *hess_d2y_i);
  _MM_DOT4_PS (da, dbPr, *hess_dxdy_r);
  _MM_DOT4_PS (da, dbPi, *hess_dxdy_i);
  _MM_DOT4_PS (da, dbPr, *hess_dydx_r);
  _MM_DOT4_PS (da, dbPi, *hess_dydx_i);
  
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  hess[0] *= dxInv*dxInv;
  hess[1] *= dxInv*dyInv;
  hess[2] *= dxInv*dyInv;
  hess[3] *= dyInv*dyInv;
#undef P
}




/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_c (UBspline_3d_c * restrict spline, 
		    double x, double y, double z,
		    complex_float* restrict val)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);

  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  // intpart = trunc (ux, uy, uz)
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  __m128i ixiyiz;
  _mm_storeu_si128 (&ixiyiz, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txtytz, txtytz);
  __m128 t3     = _mm_mul_ps (t2, txtytz);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txtytz;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, c, cPr[4], cPi[4], bcPr, bcPi,
    tmp0, tmp1, r0, r1, r2, r3, i0, i1, i2, i3;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpz,   c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0,0));  
  tmp1 = _mm_loadu_ps (P(0,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,1,0));  tmp1 = _mm_loadu_ps (P(0,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,2,0));  tmp1 = _mm_loadu_ps (P(0,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,3,0));  tmp1 = _mm_loadu_ps (P(0,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0,0));  tmp1 = _mm_loadu_ps (P(1,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,1,0));  tmp1 = _mm_loadu_ps (P(1,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,2,0));  tmp1 = _mm_loadu_ps (P(1,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,3,0));  tmp1 = _mm_loadu_ps (P(1,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0,0));  tmp1 = _mm_loadu_ps (P(2,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,1,0));  tmp1 = _mm_loadu_ps (P(2,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,2,0));  tmp1 = _mm_loadu_ps (P(2,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,3,0));  tmp1 = _mm_loadu_ps (P(2,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0,0));  tmp1 = _mm_loadu_ps (P(3,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,1,0));  tmp1 = _mm_loadu_ps (P(3,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,2,0));  tmp1 = _mm_loadu_ps (P(3,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,3,0));  tmp1 = _mm_loadu_ps (P(3,3,2));
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
eval_UBspline_3d_c_vg (UBspline_3d_c * restrict spline, 
		       double x, double y, double z,
		       complex_float* restrict val, 
		       complex_float* restrict grad)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);

  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  // intpart = trunc (ux, uy, uz)
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  __m128i ixiyiz;
  _mm_storeu_si128 (&ixiyiz, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txtytz, txtytz);
  __m128 t3     = _mm_mul_ps (t2, txtytz);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txtytz;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, c, da, db, dc,
    cPr[4], dcPr[4], bcPr, dbcPr, bdcPr,
    cPi[4], dcPi[4], bcPi, dbcPi, bdcPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpx,  da);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpy,  db);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpz,   c);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpz,  dc);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0,0));  
  tmp1 = _mm_loadu_ps (P(0,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,1,0));  tmp1 = _mm_loadu_ps (P(0,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,2,0));  tmp1 = _mm_loadu_ps (P(0,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,3,0));  tmp1 = _mm_loadu_ps (P(0,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0,0));  tmp1 = _mm_loadu_ps (P(1,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,1,0));  tmp1 = _mm_loadu_ps (P(1,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,2,0));  tmp1 = _mm_loadu_ps (P(1,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,3,0));  tmp1 = _mm_loadu_ps (P(1,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0,0));  tmp1 = _mm_loadu_ps (P(2,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,1,0));  tmp1 = _mm_loadu_ps (P(2,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,2,0));  tmp1 = _mm_loadu_ps (P(2,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,3,0));  tmp1 = _mm_loadu_ps (P(2,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0,0));  tmp1 = _mm_loadu_ps (P(3,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,1,0));  tmp1 = _mm_loadu_ps (P(3,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,2,0));  tmp1 = _mm_loadu_ps (P(3,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,3,0));  tmp1 = _mm_loadu_ps (P(3,3,2));
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


  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;

#undef P
}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_c_vgl (UBspline_3d_c * restrict spline, 
			double x, double y, double z,
			complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);

  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  // intpart = trunc (ux, uy, uz)
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  __m128i ixiyiz;
  _mm_storeu_si128 (&ixiyiz, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txtytz, txtytz);
  __m128 t3     = _mm_mul_ps (t2, txtytz);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txtytz;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cPr[4], dcPr[4], d2cPr[4], bcPr, dbcPr, bdcPr, d2bcPr, dbdcPr, bd2cPr,
    cPi[4], dcPi[4], d2cPi[4], bcPi, dbcPi, bdcPi, d2bcPi, dbdcPi, bd2cPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpz, d2c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0,0));  
  tmp1 = _mm_loadu_ps (P(0,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,1,0));  tmp1 = _mm_loadu_ps (P(0,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,2,0));  tmp1 = _mm_loadu_ps (P(0,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,3,0));  tmp1 = _mm_loadu_ps (P(0,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0,0));  tmp1 = _mm_loadu_ps (P(1,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,1,0));  tmp1 = _mm_loadu_ps (P(1,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,2,0));  tmp1 = _mm_loadu_ps (P(1,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,3,0));  tmp1 = _mm_loadu_ps (P(1,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0,0));  tmp1 = _mm_loadu_ps (P(2,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,1,0));  tmp1 = _mm_loadu_ps (P(2,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,2,0));  tmp1 = _mm_loadu_ps (P(2,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,3,0));  tmp1 = _mm_loadu_ps (P(2,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0,0));  tmp1 = _mm_loadu_ps (P(3,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,1,0));  tmp1 = _mm_loadu_ps (P(3,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,2,0));  tmp1 = _mm_loadu_ps (P(3,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,3,0));  tmp1 = _mm_loadu_ps (P(3,3,2));
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
  // Compute laplacian
  float sec_deriv[6];
  _MM_DOT4_PS (d2a, bcPr, sec_deriv[0]);
  _MM_DOT4_PS (d2a, bcPi, sec_deriv[1]);
  _MM_DOT4_PS (a, d2bcPr, sec_deriv[2]);
  _MM_DOT4_PS (a, d2bcPi, sec_deriv[3]);
  _MM_DOT4_PS (a, bd2cPr, sec_deriv[4]);
  _MM_DOT4_PS (a, bd2cPi, sec_deriv[5]);

  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  sec_deriv[0] *= dxInv*dxInv;
  sec_deriv[1] *= dxInv*dxInv;
  sec_deriv[2] *= dyInv*dyInv;
  sec_deriv[3] *= dyInv*dyInv;
  sec_deriv[4] *= dzInv*dzInv;
  sec_deriv[5] *= dzInv*dzInv;
#ifdef __cplusplus
  *lapl = std::complex<float>(sec_deriv[0] + sec_deriv[2] + sec_deriv[4],
			      sec_deriv[1] + sec_deriv[3] + sec_deriv[5]);
#else
  *lapl = (sec_deriv[0] + sec_deriv[2] + sec_deriv[4]) +
    1.0fi*(sec_deriv[1] + sec_deriv[3] + sec_deriv[5]);
#endif

#undef P
}



/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_c_vgh (UBspline_3d_c * restrict spline, 
			double x, double y, double z,
			complex_float* restrict val, 
			complex_float* restrict grad, 
			complex_float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);

  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  // intpart = trunc (ux, uy, uz)
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  __m128i ixiyiz;
  _mm_storeu_si128 (&ixiyiz, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) ((const float*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz)+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txtytz, txtytz);
  __m128 t3     = _mm_mul_ps (t2, txtytz);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txtytz;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[2], A_s[3]
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cPr[4], dcPr[4], d2cPr[4], bcPr, dbcPr, bdcPr, d2bcPr, dbdcPr, bd2cPr,
    cPi[4], dcPi[4], d2cPi[4], bcPi, dbcPi, bdcPi, d2bcPi, dbdcPi, bd2cPi,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, r0, r1, r2, r3,
    i0, i1, i2, i3;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpz, d2c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0,0));  
  tmp1 = _mm_loadu_ps (P(0,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,1,0));  tmp1 = _mm_loadu_ps (P(0,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,2,0));  tmp1 = _mm_loadu_ps (P(0,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(0,3,0));  tmp1 = _mm_loadu_ps (P(0,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[0]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[0]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0,0));  tmp1 = _mm_loadu_ps (P(1,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,1,0));  tmp1 = _mm_loadu_ps (P(1,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,2,0));  tmp1 = _mm_loadu_ps (P(1,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(1,3,0));  tmp1 = _mm_loadu_ps (P(1,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[1]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[1]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0,0));  tmp1 = _mm_loadu_ps (P(2,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,1,0));  tmp1 = _mm_loadu_ps (P(2,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,2,0));  tmp1 = _mm_loadu_ps (P(2,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(2,3,0));  tmp1 = _mm_loadu_ps (P(2,3,2));
  r3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i3   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  _MM_MATVEC4_PS (r0, r1, r2, r3,   c,   cPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3,  dc,  dcPr[2]);
  _MM_MATVEC4_PS (r0, r1, r2, r3, d2c, d2cPr[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,   c,   cPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3,  dc,  dcPi[2]);
  _MM_MATVEC4_PS (i0, i1, i2, i3, d2c, d2cPi[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0,0));  tmp1 = _mm_loadu_ps (P(3,0,2));
  r0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i0   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,1,0));  tmp1 = _mm_loadu_ps (P(3,1,2));
  r1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i1   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,2,0));  tmp1 = _mm_loadu_ps (P(3,2,2));
  r2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (2, 0, 2, 0));
  i2   = _mm_shuffle_ps (tmp0, tmp1, _MM_SHUFFLE (3, 1, 3, 1));
  tmp0 = _mm_loadu_ps (P(3,3,0));  tmp1 = _mm_loadu_ps (P(3,3,2));
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
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  hess[0] *= dxInv*dxInv;
  hess[4] *= dyInv*dyInv;
  hess[8] *= dzInv*dzInv;
  hess[1] *= dxInv*dyInv;
  hess[2] *= dxInv*dzInv;
  hess[5] *= dyInv*dzInv;
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];

#undef P
}

#undef _MM_MATVEC4_PS
#undef _MM_DOT4_PS

#endif
