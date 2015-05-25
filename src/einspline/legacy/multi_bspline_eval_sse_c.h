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

#ifndef MULTI_BSPLINE_EVAL_SSE_C_H
#define MULTI_BSPLINE_EVAL_SSE_C_H

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <math.h>
extern __m128 *restrict A_f;
extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;
#include "multi_bspline_structs.h"

#ifndef _MM_DDOT4_PD
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
  _mm_store_d (&(p), t1);                                            \
 } while (0);
#endif
#endif

/************************************************************/
/* 1D single-precision, complex evaulation functions        */
/************************************************************/
inline void
eval_multi_UBspline_1d_c (multi_UBspline_1d_c *spline,
                          double x,
                          complex_float* restrict vals)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  float tpx[4], a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
    vals[n]  = 0.0;
  for (int i=0; i<4; i++)
  {
    complex_float* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
      vals[n]  +=   a[i] * coefs[n];
  }
}



inline void
eval_multi_UBspline_1d_c_vg (multi_UBspline_1d_c *spline,
                             double x,
                             complex_float* restrict vals,
                             complex_float* restrict grads)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  float tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    complex_float* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
    {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
    }
  }
  float dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
    grads[n] *= dxInv;
}


inline void
eval_multi_UBspline_1d_c_vgl (multi_UBspline_1d_c *spline,
                              double x,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict lapl)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  float tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  complex_float* restrict coefs = spline->coefs;
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    complex_float* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
    {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
      lapl[n]  += d2a[i] * coefs[n];
    }
  }
  float dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
    grads[n] *= dxInv;
    lapl [n] *= dxInv*dxInv;
  }
}


inline void
eval_multi_UBspline_1d_c_vgh (multi_UBspline_1d_c *spline,
                              double x,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict hess)
{
  eval_multi_UBspline_1d_c_vgl (spline, x, vals, grads, hess);
}


/************************************************************/
/* 2D single-precision, complex evaulation functions        */
/************************************************************/
inline void
eval_multi_UBspline_2d_c (multi_UBspline_2d_c *spline,
                          double x, double y,
                          complex_float* restrict vals)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
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
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a4, b4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  __m128 a[4], b[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm];
  // Zero out values;
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  // Main compute loop
  __m128 ab;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],  b[j]);
      __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nm; n++)
        mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
    }
  // Now, store results back
  for (int n=0; n<N; n++)
    vals[n]      =  ((complex_float*)mvals)[n];
}


inline void
eval_multi_UBspline_2d_c_vg (multi_UBspline_2d_c *spline,
                             double x, double y,
                             complex_float* restrict vals,
                             complex_float* restrict grads)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
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
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a4, b4, da4, db4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  __m128 a[4], b[4], da[4], db[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[2*Nm];
  // Zero out values;
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  // Main compute loop
  __m128 ab, dab[2];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],  b[j]);
      dab[0]  = _mm_mul_ps ( da[i],  b[j]);
      dab[1]  = _mm_mul_ps (  a[i], db[j]);
      __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nm; n++)
      {
        mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
        mgrad[2*n+0] = _mm_add_ps (mgrad[2*n+0], _mm_mul_ps( dab[0], coefs[n]));
        mgrad[2*n+1] = _mm_add_ps (mgrad[2*n+1], _mm_mul_ps( dab[1], coefs[n]));
      }
    }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      =  ((complex_float*)mvals)[n];
    grads[2*n+0] =  ((complex_float*)mgrad)[nd2*4 + 2*0 + nm2] * dxInv;
    grads[2*n+1] =  ((complex_float*)mgrad)[nd2*4 + 2*1 + nm2] * dyInv;
  }
}

inline void
eval_multi_UBspline_2d_c_vgl (multi_UBspline_2d_c *spline,
                              double x, double y,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);
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
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a4, b4, da4, db4, d2a4, d2b4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b4);
  __m128 a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2a4, d2a4);
  d2a[0]=_mm_unpacklo_ps(tmp, tmp);
  d2a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2a4, d2a4);
  d2a[2]=_mm_unpacklo_ps(tmp, tmp);
  d2a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2b4, d2b4);
  d2b[0]=_mm_unpacklo_ps(tmp, tmp);
  d2b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2b4, d2b4);
  d2b[2]=_mm_unpacklo_ps(tmp, tmp);
  d2b[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[2*Nm], mlapl[2*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mlapl[n] = _mm_setzero_ps();
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  // Main compute loop
  __m128 ab, dab[2], d2ab[2];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],  b[j]);
      dab[0]  = _mm_mul_ps ( da[i],  b[j]);
      dab[1]  = _mm_mul_ps (  a[i], db[j]);
      d2ab[0] = _mm_mul_ps (d2a[i],  b[j]);
      d2ab[1] = _mm_mul_ps (  a[i],d2b[j]);
      __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nm; n++)
      {
        mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
        mgrad[2*n+0] = _mm_add_ps (mgrad[2*n+0], _mm_mul_ps( dab[0], coefs[n]));
        mgrad[2*n+1] = _mm_add_ps (mgrad[2*n+1], _mm_mul_ps( dab[1], coefs[n]));
        mlapl[2*n+0] = _mm_add_ps (mlapl[2*n+0], _mm_mul_ps(d2ab[0], coefs[n]));
        mlapl[2*n+1] = _mm_add_ps (mlapl[2*n+1], _mm_mul_ps(d2ab[1], coefs[n]));
      }
    }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      =  ((complex_float*)mvals)[n];
    grads[2*n+0] =  ((complex_float*)mgrad)[nd2*4 + 2*0 + nm2] * dxInv;
    grads[2*n+1] =  ((complex_float*)mgrad)[nd2*4 + 2*1 + nm2] * dyInv;
    lapl [n]     = (((complex_float*)mlapl)[nd2*4 + 2*0 + nm2] * dxInv*dxInv +
                    ((complex_float*)mlapl)[nd2*4 + 2*1 + nm2] * dyInv*dyInv);
  }
}

inline void
eval_multi_UBspline_2d_c_vgh (multi_UBspline_2d_c *spline,
                              double x, double y,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);
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
//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]
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
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  __m128 tpx = _mm_set_ps (tx*tx*tx, tx*tx, tx, 1.0);
  __m128 tpy = _mm_set_ps (ty*ty*ty, ty*ty, ty, 1.0);
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3]
  __m128 a4, b4, da4, db4, d2a4, d2b4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b4);
  __m128 a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2a4, d2a4);
  d2a[0]=_mm_unpacklo_ps(tmp, tmp);
  d2a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2a4, d2a4);
  d2a[2]=_mm_unpacklo_ps(tmp, tmp);
  d2a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2b4, d2b4);
  d2b[0]=_mm_unpacklo_ps(tmp, tmp);
  d2b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2b4, d2b4);
  d2b[2]=_mm_unpacklo_ps(tmp, tmp);
  d2b[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[2*Nm], mhess[3*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mhess[n] = _mm_setzero_ps();
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  // Main compute loop
  __m128 ab, dab[2], d2ab[3];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],  b[j]);
      dab[0]  = _mm_mul_ps ( da[i],  b[j]);
      dab[1]  = _mm_mul_ps (  a[i], db[j]);
      d2ab[0] = _mm_mul_ps (d2a[i],  b[j]);
      d2ab[1] = _mm_mul_ps ( da[i], db[j]);
      d2ab[2] = _mm_mul_ps (  a[i],d2b[j]);
      __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nm; n++)
      {
        mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
        mgrad[2*n+0] = _mm_add_ps (mgrad[2*n+0], _mm_mul_ps( dab[0], coefs[n]));
        mgrad[2*n+1] = _mm_add_ps (mgrad[2*n+1], _mm_mul_ps( dab[1], coefs[n]));
        mhess[3*n+0] = _mm_add_ps (mhess[3*n+0], _mm_mul_ps(d2ab[0], coefs[n]));
        mhess[3*n+1] = _mm_add_ps (mhess[3*n+1], _mm_mul_ps(d2ab[1], coefs[n]));
        mhess[3*n+2] = _mm_add_ps (mhess[3*n+2], _mm_mul_ps(d2ab[2], coefs[n]));
      }
    }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      = ((complex_float*)mvals)[n];
    grads[2*n+0] = ((complex_float*)mgrad)[nd2*4 + 2*0 + nm2] * dxInv;
    grads[2*n+1] = ((complex_float*)mgrad)[nd2*4 + 2*1 + nm2] * dyInv;
    hess [4*n+0]               = ((complex_float*)mhess)[nd2*6 + 2*0 + nm2] * dxInv*dxInv;
    hess [4*n+1] = hess[4*n+2] = ((complex_float*)mhess)[nd2*6 + 2*1 + nm2] * dxInv*dyInv;
    hess [4*n+3]               = ((complex_float*)mhess)[nd2*6 + 2*2 + nm2] * dyInv*dyInv;
  }
}


/************************************************************/
/* 3D single-precision, complex evaulation functions        */
/************************************************************/
inline void
eval_multi_UBspline_3d_c (multi_UBspline_3d_c *spline,
                          double x, double y, double z,
                          complex_float* restrict vals)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
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
  __m128 a4, b4, c4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
  __m128 a[4], b[4], c[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack c values
  tmp=_mm_unpacklo_ps(  c4,   c4);
  c[0]=_mm_unpacklo_ps(tmp, tmp);
  c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  c4,   c4);
  c[2]=_mm_unpacklo_ps(tmp, tmp);
  c[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm];
  // Zero out values;
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  // Main compute loop
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  __m128 abc;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
          mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  abc   , coefs[n]));
      }
  // Now, store results back
  for (int n=0; n<N; n++)
    vals[n]      = ((complex_float*)mvals)[n];
}



inline void
eval_multi_UBspline_3d_c_vg (multi_UBspline_3d_c *spline,
                             double x, double y, double z,
                             complex_float* restrict vals,
                             complex_float* restrict grads)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc4);
  __m128 a[4], b[4], c[4], da[4], db[4], dc[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack c values
  tmp=_mm_unpacklo_ps(  c4,   c4);
  c[0]=_mm_unpacklo_ps(tmp, tmp);
  c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  c4,   c4);
  c[2]=_mm_unpacklo_ps(tmp, tmp);
  c[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( dc4,  dc4);
  dc[0]=_mm_unpacklo_ps(tmp, tmp);
  dc[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( dc4,  dc4);
  dc[2]=_mm_unpacklo_ps(tmp, tmp);
  dc[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[3*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 abc, dabc[3];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[0]  = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[1]  = _mm_mul_ps (  a[i], _mm_mul_ps( db[j],  c[k]));
        dabc[2]  = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  dc[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
        {
          mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  abc   , coefs[n]));
          mgrad[3*n+0] = _mm_add_ps (mgrad[3*n+0], _mm_mul_ps( dabc[0], coefs[n]));
          mgrad[3*n+1] = _mm_add_ps (mgrad[3*n+1], _mm_mul_ps( dabc[1], coefs[n]));
          mgrad[3*n+2] = _mm_add_ps (mgrad[3*n+2], _mm_mul_ps( dabc[2], coefs[n]));
        }
      }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      = ((complex_float*)mvals)[n];
    grads[3*n+0] = ((complex_float*)mgrad)[nd2*6 + 2*0 + nm2] * dxInv;
    grads[3*n+1] = ((complex_float*)mgrad)[nd2*6 + 2*1 + nm2] * dyInv;
    grads[3*n+2] = ((complex_float*)mgrad)[nd2*6 + 2*2 + nm2] * dzInv;
  }
}



inline void
eval_multi_UBspline_3d_c_vgl (multi_UBspline_3d_c *spline,
                              double x, double y, double z,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b4);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpz, d2c4);
  __m128 a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2a4, d2a4);
  d2a[0]=_mm_unpacklo_ps(tmp, tmp);
  d2a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2a4, d2a4);
  d2a[2]=_mm_unpacklo_ps(tmp, tmp);
  d2a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2b4, d2b4);
  d2b[0]=_mm_unpacklo_ps(tmp, tmp);
  d2b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2b4, d2b4);
  d2b[2]=_mm_unpacklo_ps(tmp, tmp);
  d2b[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack c values
  tmp=_mm_unpacklo_ps(  c4,   c4);
  c[0]=_mm_unpacklo_ps(tmp, tmp);
  c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  c4,   c4);
  c[2]=_mm_unpacklo_ps(tmp, tmp);
  c[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( dc4,  dc4);
  dc[0]=_mm_unpacklo_ps(tmp, tmp);
  dc[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( dc4,  dc4);
  dc[2]=_mm_unpacklo_ps(tmp, tmp);
  dc[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2c4, d2c4);
  d2c[0]=_mm_unpacklo_ps(tmp, tmp);
  d2c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2c4, d2c4);
  d2c[2]=_mm_unpacklo_ps(tmp, tmp);
  d2c[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[3*Nm], mlapl[3*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mlapl[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 abc, dabc[3], d2abc[3];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[0]  = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[1]  = _mm_mul_ps (  a[i], _mm_mul_ps( db[j],  c[k]));
        dabc[2]  = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  dc[k]));
        d2abc[0] = _mm_mul_ps (d2a[i], _mm_mul_ps(  b[j],  c[k]));
        d2abc[1] = _mm_mul_ps (  a[i], _mm_mul_ps(d2b[j],  c[k]));
        d2abc[2] = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j], d2c[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
        {
          mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  abc   , coefs[n]));
          mgrad[3*n+0] = _mm_add_ps (mgrad[3*n+0], _mm_mul_ps( dabc[0], coefs[n]));
          mgrad[3*n+1] = _mm_add_ps (mgrad[3*n+1], _mm_mul_ps( dabc[1], coefs[n]));
          mgrad[3*n+2] = _mm_add_ps (mgrad[3*n+2], _mm_mul_ps( dabc[2], coefs[n]));
          mlapl[3*n+0] = _mm_add_ps (mlapl[3*n+0], _mm_mul_ps(d2abc[0], coefs[n]));
          mlapl[3*n+1] = _mm_add_ps (mlapl[3*n+1], _mm_mul_ps(d2abc[1], coefs[n]));
          mlapl[3*n+2] = _mm_add_ps (mlapl[3*n+2], _mm_mul_ps(d2abc[2], coefs[n]));
        }
      }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      = ((complex_float*)mvals)[n];
    grads[3*n+0] = ((complex_float*)mgrad)[nd2*6 + 2*0 + nm2] * dxInv;
    grads[3*n+1] = ((complex_float*)mgrad)[nd2*6 + 2*1 + nm2] * dyInv;
    grads[3*n+2] = ((complex_float*)mgrad)[nd2*6 + 2*2 + nm2] * dzInv;
    lapl [n] = (((complex_float*)mlapl)[nd2*6 + 2*0 + nm2] * dxInv*dxInv +
                ((complex_float*)mlapl)[nd2*6 + 2*1 + nm2] * dyInv*dyInv +
                ((complex_float*)mlapl)[nd2*6 + 2*2 + nm2] * dzInv*dzInv);
  }
}


inline void
eval_multi_UBspline_3d_c_vgh (multi_UBspline_3d_c *spline,
                              double x, double y, double z,
                              complex_float* restrict vals,
                              complex_float* restrict grads,
                              complex_float* restrict hess)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[10],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[11],_MM_HINT_T0);
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b4);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpz, d2c4);
  __m128 a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  __m128 tmp;
  // Unpack a values
  tmp=_mm_unpacklo_ps(  a4,   a4);
  a[0]=_mm_unpacklo_ps(tmp, tmp);
  a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  a4,   a4);
  a[2]=_mm_unpacklo_ps(tmp, tmp);
  a[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( da4,  da4);
  da[0]=_mm_unpacklo_ps(tmp, tmp);
  da[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( da4,  da4);
  da[2]=_mm_unpacklo_ps(tmp, tmp);
  da[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2a4, d2a4);
  d2a[0]=_mm_unpacklo_ps(tmp, tmp);
  d2a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2a4, d2a4);
  d2a[2]=_mm_unpacklo_ps(tmp, tmp);
  d2a[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack b values
  tmp=_mm_unpacklo_ps(  b4,   b4);
  b[0]=_mm_unpacklo_ps(tmp, tmp);
  b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  b4,   b4);
  b[2]=_mm_unpacklo_ps(tmp, tmp);
  b[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( db4,  db4);
  db[0]=_mm_unpacklo_ps(tmp, tmp);
  db[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( db4,  db4);
  db[2]=_mm_unpacklo_ps(tmp, tmp);
  db[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2b4, d2b4);
  d2b[0]=_mm_unpacklo_ps(tmp, tmp);
  d2b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2b4, d2b4);
  d2b[2]=_mm_unpacklo_ps(tmp, tmp);
  d2b[3]=_mm_unpackhi_ps(tmp, tmp);
  // Unpack c values
  tmp=_mm_unpacklo_ps(  c4,   c4);
  c[0]=_mm_unpacklo_ps(tmp, tmp);
  c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(  c4,   c4);
  c[2]=_mm_unpacklo_ps(tmp, tmp);
  c[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps( dc4,  dc4);
  dc[0]=_mm_unpacklo_ps(tmp, tmp);
  dc[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps( dc4,  dc4);
  dc[2]=_mm_unpacklo_ps(tmp, tmp);
  dc[3]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpacklo_ps(d2c4, d2c4);
  d2c[0]=_mm_unpacklo_ps(tmp, tmp);
  d2c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d2c4, d2c4);
  d2c[2]=_mm_unpacklo_ps(tmp, tmp);
  d2c[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+1)/2;
  __m128 mvals[Nm], mgrad[3*Nm], mhess[6*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<6*Nm; n++)
    mhess[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 abc, dabc[3], d2abc[6];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[0]  = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],  c[k]));
        dabc[1]  = _mm_mul_ps (  a[i], _mm_mul_ps( db[j],  c[k]));
        dabc[2]  = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  dc[k]));
        d2abc[0] = _mm_mul_ps (d2a[i], _mm_mul_ps(  b[j],  c[k]));
        d2abc[1] = _mm_mul_ps ( da[i], _mm_mul_ps( db[j],  c[k]));
        d2abc[2] = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],  dc[k]));
        d2abc[3] = _mm_mul_ps (  a[i], _mm_mul_ps(d2b[j],  c[k]));
        d2abc[4] = _mm_mul_ps (  a[i], _mm_mul_ps( db[j], dc[k]));
        d2abc[5] = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j], d2c[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
        {
          mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  abc   , coefs[n]));
          mgrad[3*n+0] = _mm_add_ps (mgrad[3*n+0], _mm_mul_ps( dabc[0], coefs[n]));
          mgrad[3*n+1] = _mm_add_ps (mgrad[3*n+1], _mm_mul_ps( dabc[1], coefs[n]));
          mgrad[3*n+2] = _mm_add_ps (mgrad[3*n+2], _mm_mul_ps( dabc[2], coefs[n]));
          mhess[6*n+0] = _mm_add_ps (mhess[6*n+0], _mm_mul_ps(d2abc[0], coefs[n]));
          mhess[6*n+1] = _mm_add_ps (mhess[6*n+1], _mm_mul_ps(d2abc[1], coefs[n]));
          mhess[6*n+2] = _mm_add_ps (mhess[6*n+2], _mm_mul_ps(d2abc[2], coefs[n]));
          mhess[6*n+3] = _mm_add_ps (mhess[6*n+3], _mm_mul_ps(d2abc[3], coefs[n]));
          mhess[6*n+4] = _mm_add_ps (mhess[6*n+4], _mm_mul_ps(d2abc[4], coefs[n]));
          mhess[6*n+5] = _mm_add_ps (mhess[6*n+5], _mm_mul_ps(d2abc[5], coefs[n]));
        }
      }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd2 = n>>1;
    int nm2 = n & 1;
    vals[n]      = ((complex_float*)mvals)[n];
    grads[3*n+0] = ((complex_float*)mgrad)[nd2*6 + 2*0 + nm2] * dxInv;
    grads[3*n+1] = ((complex_float*)mgrad)[nd2*6 + 2*1 + nm2] * dyInv;
    grads[3*n+2] = ((complex_float*)mgrad)[nd2*6 + 2*2 + nm2] * dzInv;
    hess [9*n+0]               = ((complex_float*)mhess)[nd2*12 + 2*0 + nm2] * dxInv*dxInv;
    hess [9*n+1] = hess[9*n+3] = ((complex_float*)mhess)[nd2*12 + 2*1 + nm2] * dxInv*dyInv;
    hess [9*n+2] = hess[9*n+6] = ((complex_float*)mhess)[nd2*12 + 2*2 + nm2] * dxInv*dzInv;
    hess [9*n+4]               = ((complex_float*)mhess)[nd2*12 + 2*3 + nm2] * dyInv*dyInv;
    hess [9*n+5] = hess[9*n+7] = ((complex_float*)mhess)[nd2*12 + 2*4 + nm2] * dyInv*dzInv;
    hess [9*n+8]               = ((complex_float*)mhess)[nd2*12 + 2*5 + nm2] * dzInv*dzInv;
  }
}

#endif
