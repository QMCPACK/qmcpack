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

#ifndef BSPLINE_EVAL_SSE_S_H
#define BSPLINE_EVAL_SSE_S_H

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#include <stdio.h>
#include <math.h>

// extern __m128   A0,   A1,   A2,   A3;
// extern __m128  dA0,  dA1,  dA2,  dA3;
// extern __m128 d2A0, d2A1, d2A2, d2A3;
extern __m128* restrict A_s;

extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;

/// SSE3 add "horizontal add" instructions, which makes things
/// simpler and faster
#ifdef HAVE_SSE9
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 r0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
  __m128 r1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
  r = _mm_hadd_ps (r0, r1);                                         \
 } while (0);
#define _MM_DOT4_PS(A, B, _p)                                       \
do {                                                                \
  __m128 t  = _mm_mul_ps (A, B);                                    \
  __m128 t1 = _mm_hadd_ps (t,t);                                    \
  __m128 r  = _mm_hadd_ps (t1, t1);                                 \
  _mm_store_ss (&(_p), r);                                          \
} while(0);
#else
// Use plain-old SSE instructions
#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  __m128 _r0 = _mm_mul_ps (M0, v);                                  \
  __m128 _r1 = _mm_mul_ps (M1, v);				    \
  __m128 _r2 = _mm_mul_ps (M2, v);                                  \
  __m128 _r3 = _mm_mul_ps (M3, v);				    \
  _MM_TRANSPOSE4_PS (_r0, _r1, _r2, _r3);                           \
  r = _mm_add_ps (_mm_add_ps (_r0, _r1), _mm_add_ps (_r2, _r3));    \
 } while (0);
#define _MM_DOT4_PS(A, B, p)                                        \
do {                                                                \
  __m128 t    = _mm_mul_ps (A, B);                                  \
  __m128 alo  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(0,1,0,1));	    \
  __m128 ahi  = _mm_shuffle_ps (t, t, _MM_SHUFFLE(2,3,2,3));	    \
  __m128 _a    = _mm_add_ps (alo, ahi);                              \
  __m128 rlo  = _mm_shuffle_ps (_a, _a, _MM_SHUFFLE(0,0,0,0));	    \
  __m128 rhi  = _mm_shuffle_ps (_a, _a, _MM_SHUFFLE(1,1,1,1));	    \
  __m128 _r   = _mm_add_ps (rlo, rhi);                              \
  _mm_store_ss (&(p), _r);                                          \
} while(0);
#endif


/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_s (UBspline_1d_s * restrict spline, 
		    double x, float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x, 
		     float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x, 
			float* restrict val, float* restrict grad,
			float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float* restrict coefs = spline->coefs;
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;

  // It turns out that std version is faster than SSE in 1D
//   __m128 tp = _mm_set_ps (t*t*t, t*t, t, 1.0);
//   __m128 a, da, d2a;
//   _MM_MATVEC4_PS (  A_s[0],   A_s[1],   A_s[ 2],   A_s[ 3], tp,   a);
//   _MM_MATVEC4_PS ( A_s[4],  A_s[5],  A_s[6],  A_s[7], tp,  da);
//   _MM_MATVEC4_PS (A_s[8], A_s[9], A_s[10], A_s[11], tp, d2a);
//   __m128 cf  = _mm_loadu_ps (&(coefs[i]));

//   _MM_DOT4_PS (  a, cf, *val);
//   _MM_DOT4_PS ( da, cf, *grad);
//   _MM_DOT4_PS (d2a, cf, *lapl);


  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
     coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
     coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
     coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

inline void
eval_UBspline_1d_s_vgh (UBspline_1d_s * restrict spline, double x, 
			float* restrict val, float* restrict grad,
			float* restrict hess)
{
  eval_UBspline_1d_s_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
		    double x, double y, float* restrict val)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
  /// SSE mesh point determination
  __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
  __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
  xy = _mm_sub_ps (xy, x0y0);
  // ux = (x - x0)/delta_x and same for y
  __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
  // intpart = trunc (ux, uy)
  __m128i intpart  = _mm_cvttps_epi32(uxuy);
  __m128i ixiy;
  _mm_storeu_si128 (&ixiy, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiy)[3];
  int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txty   = _mm_sub_ps (uxuy, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txty, txty);
  __m128 t3     = _mm_mul_ps (t2, txty);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txty;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, bP, tmp0, tmp1, tmp2, tmp3;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0));
  tmp1 = _mm_loadu_ps (P(1));
  tmp2 = _mm_loadu_ps (P(2));
  tmp3 = _mm_loadu_ps (P(3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
  // Compute value
  _MM_DOT4_PS (a, bP, *val);
#undef P
}


/* Value and gradient */
inline void
eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
		       double x, double y, 
		       float* restrict val, float* restrict grad)
{
  _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[7],_MM_HINT_T0);
  /// SSE mesh point determination
  __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
  __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
  xy = _mm_sub_ps (xy, x0y0);
  // ux = (x - x0)/delta_x and same for y
  __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
  // intpart = trunc (ux, uy)
  __m128i intpart  = _mm_cvttps_epi32(uxuy);
  __m128i ixiy;
  _mm_storeu_si128 (&ixiy, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiy)[3];
  int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txty   = _mm_sub_ps (uxuy, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txty, txty);
  __m128 t3     = _mm_mul_ps (t2, txty);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txty;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, da, db, bP, dbP, tmp0, tmp1, tmp2, tmp3;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpx,  da);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpy,  db);
  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0));
  tmp1 = _mm_loadu_ps (P(1));
  tmp2 = _mm_loadu_ps (P(2));
  tmp3 = _mm_loadu_ps (P(3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  db,  dbP);
  // Compute value
  _MM_DOT4_PS (a, bP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bP, grad[0]);
  _MM_DOT4_PS (a, dbP, grad[1]);
  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);  
  /// SSE mesh point determination
  __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
  __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
  xy = _mm_sub_ps (xy, x0y0);
  // ux = (x - x0)/delta_x and same for y
  __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
  // intpart = trunc (ux, uy)
  __m128i intpart  = _mm_cvttps_epi32(uxuy);
  __m128i ixiy;
  _mm_storeu_si128 (&ixiy, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiy)[3];
  int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txty   = _mm_sub_ps (uxuy, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txty, txty);
  __m128 t3     = _mm_mul_ps (t2, txty);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txty;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, da, db, d2a, d2b, bP, dbP, d2bP, 
    tmp0, tmp1, tmp2, tmp3;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[8], A_s[9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[8], A_s[9], A_s[10], A_s[11], tpy, d2b);
  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0));
  tmp1 = _mm_loadu_ps (P(1));
  tmp2 = _mm_loadu_ps (P(2));
  tmp3 = _mm_loadu_ps (P(3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  db,  dbP);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2b, d2bP);
  // Compute value
  _MM_DOT4_PS (a, bP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bP, grad[0]);
  _MM_DOT4_PS (a, dbP, grad[1]);
  float sec_derivs[2];
  // Compute laplacian
  _MM_DOT4_PS (d2a, bP, sec_derivs[0]);
  _MM_DOT4_PS (a, d2bP, sec_derivs[1]);

  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  sec_derivs[0] *= dxInv*dxInv;
  sec_derivs[1] *= dyInv*dyInv;
  // Copy hessian elements into lower half of 2x2 matrix
  *lapl = sec_derivs[0] + sec_derivs[1];
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline, 
			double x, double y, float* restrict val, 
			float* restrict grad, float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);  
  /// SSE mesh point determination
  __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
  __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
  xy = _mm_sub_ps (xy, x0y0);
  // ux = (x - x0)/delta_x and same for y
  __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
  // intpart = trunc (ux, uy)
  __m128i intpart  = _mm_cvttps_epi32(uxuy);
  __m128i ixiy;
  _mm_storeu_si128 (&ixiy, intpart);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiy)[3];
  int iy = ((int *)&ixiy)[2];

  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no j value is needed.
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  __m128 ipart  = _mm_cvtepi32_ps (intpart);
  __m128 txty   = _mm_sub_ps (uxuy, ipart);
  __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
  __m128 t2     = _mm_mul_ps (txty, txty);
  __m128 t3     = _mm_mul_ps (t2, txty);
  __m128 tpx    = t3;
  __m128 tpy    = t2;
  __m128 tpz    = txty;
  __m128 zero   = one;
  _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, da, db, d2a, d2b, bP, dbP, d2bP, 
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[ 2], A_s[ 3], tpx,   a);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[ 6], A_s[ 7], tpx,  da);
  _MM_MATVEC4_PS (A_s[8], A_s[9], A_s[10], A_s[11], tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[ 2], A_s[ 3], tpy,   b);
  _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[ 6], A_s[ 7], tpy,  db);
  _MM_MATVEC4_PS (A_s[8], A_s[9], A_s[10], A_s[11], tpy, d2b);
  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0));
  tmp1 = _mm_loadu_ps (P(1));
  tmp2 = _mm_loadu_ps (P(2));
  tmp3 = _mm_loadu_ps (P(3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  db,  dbP);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2b, d2bP);
  // Compute value
  _MM_DOT4_PS (a, bP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bP, grad[0]);
  _MM_DOT4_PS (a, dbP, grad[1]);
  // Compute hessian
  _MM_DOT4_PS (d2a, bP, hess[0]);
  _MM_DOT4_PS (a, d2bP, hess[3]);
  _MM_DOT4_PS (da, dbP, hess[1]);

  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  hess[0] *= dxInv*dxInv;
  hess[3] *= dyInv*dyInv;
  hess[1] *= dxInv*dyInv;
  // Copy hessian elements into lower half of 2x2 matrix
  hess[2] = hess[1];
#undef P
}




/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
		    double x, double y, double z,
		    float* restrict val)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);

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
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

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
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, c, cP[4],bcP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpz,   c);

  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);

#undef P
}

/* Value and gradient */
inline void
eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);

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
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

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
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, c, da, db, dc,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
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
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  
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
eval_UBspline_3d_s_vgl (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);  

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
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

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
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
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
  tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));
  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));
  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3], d2b, d2bcP);
  _MM_MATVEC4_PS (d2cP[0], d2cP[1], d2cP[2], d2cP[3],   b, bd2cP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],  db, dbdcP);

  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  // Compute laplacian
  float lx, ly, lz;
  _MM_DOT4_PS (d2a, bcP, lx);
  _MM_DOT4_PS (a, d2bcP, ly);
  _MM_DOT4_PS (a, bd2cP, lz);
  
  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  lx *= dxInv*dxInv;
  ly *= dyInv*dyInv;
  lz *= dzInv*dzInv;
  *lapl = lx + ly + lz;	       
#undef P
}


/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline, 
			double x, double y, double z,
			float* restrict val, float* restrict grad, 
			float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);  

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
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

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
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a, b, c, da, db, dc, d2a, d2b, d2c,
    cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

  
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
  tmp0 = _mm_loadu_ps (P(0,0));
  tmp1 = _mm_loadu_ps (P(0,1));
  tmp2 = _mm_loadu_ps (P(0,2));
  tmp3 = _mm_loadu_ps (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[0]);
  // 2nd quarter
  tmp0 = _mm_loadu_ps (P(1,0));
  tmp1 = _mm_loadu_ps (P(1,1));
  tmp2 = _mm_loadu_ps (P(1,2));
  tmp3 = _mm_loadu_ps (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[1]);
  // 3rd quarter
  tmp0 = _mm_loadu_ps (P(2,0));
  tmp1 = _mm_loadu_ps (P(2,1));
  tmp2 = _mm_loadu_ps (P(2,2));
  tmp3 = _mm_loadu_ps (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[2]);
  // 4th quarter
  tmp0 = _mm_loadu_ps (P(3,0));
  tmp1 = _mm_loadu_ps (P(3,1));
  tmp2 = _mm_loadu_ps (P(3,2));
  tmp3 = _mm_loadu_ps (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[3]);
  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3], d2b, d2bcP);
  _MM_MATVEC4_PS (d2cP[0], d2cP[1], d2cP[2], d2cP[3],   b, bd2cP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],  db, dbdcP);
  // Compute value
  _MM_DOT4_PS (a, bcP, *val);
  // Compute gradient
  _MM_DOT4_PS (da, bcP, grad[0]);
  _MM_DOT4_PS (a, dbcP, grad[1]);
  _MM_DOT4_PS (a, bdcP, grad[2]);
  // Compute hessian
  _MM_DOT4_PS (d2a, bcP, hess[0]);
  _MM_DOT4_PS (a, d2bcP, hess[4]);
  _MM_DOT4_PS (a, bd2cP, hess[8]);
  _MM_DOT4_PS (da, dbcP, hess[1]);
  _MM_DOT4_PS (da, bdcP, hess[2]);
  _MM_DOT4_PS (a, dbdcP, hess[5]);

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
