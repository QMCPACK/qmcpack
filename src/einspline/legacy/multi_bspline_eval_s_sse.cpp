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

#include <config.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#ifdef HAVE_SSE41
#include <smmintrin.h>
#endif
#include <stdio.h>
#include <math.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"
#include "multi_bspline_eval_s.h"

using std::min;
using std::max;

extern __m128 *restrict A_s;
extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;
extern const float* restrict d3Af;


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

#if !defined(HAVE_SSE41)
inline __m128i _mm_min_epi32(__m128i a, __m128i b)
{
  __m128i mask  = _mm_cmplt_epi32(a, b);
  a = _mm_and_si128(a, mask);
  b = _mm_andnot_si128(mask, b);
  a = _mm_or_si128(a, b);
  return a;
}

inline __m128i _mm_max_epi32(__m128i a, __m128i b)
{
  __m128i mask  = _mm_cmpgt_epi32(a, b);
  a = _mm_and_si128(a, mask);
  b = _mm_andnot_si128(mask, b);
  a = _mm_or_si128(a, b);
  return a;
}

#endif

/************************************************************/
/* 1D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_s (const multi_UBspline_1d_s *spline,
                          float x,
                          float* restrict vals)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  float tpx[4], a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
    vals[n]  = 0.0;
  for (int i=0; i<4; i++)
  {
    float* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
      vals[n]  +=   a[i] * coefs[n];
  }
}



void
eval_multi_UBspline_1d_s_vg (const multi_UBspline_1d_s *spline,
                             float x,
                             float* restrict vals,
                             float* restrict grads)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  float tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;
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
    float* restrict coefs = spline->coefs + ((ix+i)*xs);
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


void
eval_multi_UBspline_1d_s_vgl (const multi_UBspline_1d_s *spline,
                              float x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  float tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;
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
    float* restrict coefs = spline->coefs + ((ix+i)*xs);
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


void
eval_multi_UBspline_1d_s_vgh (const multi_UBspline_1d_s *spline,
                              float x,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess)
{
  eval_multi_UBspline_1d_s_vgl (spline, x, vals, grads, hess);
}

/************************************************************/
/* 2D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_s(const multi_UBspline_2d_s *spline,
                         double x, double y,
                         float* restrict vals)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
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
  int Nm = (N+3)/4;
  __m128 mvals[Nm];
  // Zero out values;
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 ab, dab[2];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],   b[j]);
      __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nm; n++)
        mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
    }
  // Now, store results back
  for (int n=0; n<N; n++)
    vals[n] = ((float*)mvals)[n];
}


void
eval_multi_UBspline_2d_s_vg (const multi_UBspline_2d_s *spline,
                             double x, double y,
                             float* restrict vals,
                             float* restrict grads)
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4;
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
  int Nm = (N+3)/4;
  __m128 mvals[Nm], mgrad[2*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 ab, dab[2];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],   b[j]);
      dab[0]  = _mm_mul_ps ( da[i],   b[j]);
      dab[1]  = _mm_mul_ps (  a[i],  db[j]);
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[2*n+0] = ((float*)mgrad)[nd4*8 + 4*0 + nm4] * dxInv;
    grads[2*n+1] = ((float*)mgrad)[nd4*8 + 4*1 + nm4] * dyInv;
  }
}



void
eval_multi_UBspline_2d_s_vgl (const multi_UBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl)
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4;
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
  int Nm = (N+3)/4;
  __m128 mvals[Nm], mgrad[2*Nm], mlapl[2*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mlapl[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 ab, dab[2], d2ab[2];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],   b[j]);
      dab[0]  = _mm_mul_ps ( da[i],   b[j]);
      dab[1]  = _mm_mul_ps (  a[i],  db[j]);
      d2ab[0] = _mm_mul_ps (d2a[i],   b[j]);
      d2ab[1] = _mm_mul_ps (  a[i], d2b[j]);
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[2*n+0] = ((float*)mgrad)[nd4*8 + 4*0 + nm4] * dxInv;
    grads[2*n+1] = ((float*)mgrad)[nd4*8 + 4*1 + nm4] * dyInv;
    lapl [n]     = (((float*)mlapl)[nd4*8 + 4*0 + nm4] * dxInv*dxInv +
                    ((float*)mlapl)[nd4*8 + 4*1 + nm4] * dyInv*dyInv);
  }
}



void
eval_multi_UBspline_2d_s_vgh (const multi_UBspline_2d_s *spline,
                              double x, double y,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess)
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
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4;
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
  int Nm = (N+3)/4;
  __m128 mvals[Nm], mgrad[2*Nm], mhess[3*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<2*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mhess[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 ab, dab[2], d2ab[3];
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      ab      = _mm_mul_ps (  a[i],   b[j]);
      dab[0]  = _mm_mul_ps ( da[i],   b[j]);
      dab[1]  = _mm_mul_ps (  a[i],  db[j]);
      d2ab[0] = _mm_mul_ps (d2a[i],   b[j]);
      d2ab[1] = _mm_mul_ps ( da[i],  db[j]);
      d2ab[2] = _mm_mul_ps (  a[i], d2b[j]);
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[2*n+0] = ((float*)mgrad)[nd4*8 + 4*0 + nm4] * dxInv;
    grads[2*n+1] = ((float*)mgrad)[nd4*8 + 4*1 + nm4] * dyInv;
    hess [4*n+0]               = ((float*)mhess)[nd4*12 + 4*0 + nm4] * dxInv*dxInv;
    hess [4*n+1] = hess[4*n+2] = ((float*)mhess)[nd4*12 + 4*1 + nm4] * dxInv*dyInv;
    hess [4*n+3]               = ((float*)mhess)[nd4*12 + 4*2 + nm4] * dyInv*dyInv;
  }
}


/************************************************************/
/* 3D single-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_3d_s (const multi_UBspline_3d_s *spline,
                          float x, float y, float z,
                          float* restrict vals)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
#if defined(HAVE_SSE41)
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
  //__m128i ixiyiz;
  //_mm_storeu_si128 (&ixiyiz, intpart);
  __m128i gmin = _mm_set_epi32(0,0,0,0);
  __m128i gmax = _mm_set_epi32(spline->x_grid.num-1, spline->y_grid.num-1,spline->z_grid.num-1,0);
  __m128i ixiyiz=_mm_min_epi32(_mm_max_epi32(intpart,gmin),gmax);
  _mm_storeu_si128 (&intpart, ixiyiz);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];
#else
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);
  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);
  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  __m128  uxuyuz= _mm_set_ps (ux, uy, uz, 0.0);
  __m128i intpart=_mm_set_epi32(ix,iy,iz,0);
#endif
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
  int Nm = (N+3)/4;
  __m128 mvals[Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 abc;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
          mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
      }
  // Now, store results back
  for (int n=0; n<N; n++)
    vals[n] = ((float*)mvals)[n];
}



void
eval_multi_UBspline_3d_s_vg (const multi_UBspline_3d_s *spline,
                             float x, float y, float z,
                             float* restrict vals,
                             float* restrict grads)
{
  _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
#if defined(HAVE_SSE41)
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
  //__m128i ixiyiz;
  //_mm_storeu_si128 (&ixiyiz, intpart);
  __m128i gmin = _mm_set_epi32(0,0,0,0);
  __m128i gmax = _mm_set_epi32(spline->x_grid.num-1, spline->y_grid.num-1,spline->z_grid.num-1,0);
  __m128i ixiyiz=_mm_min_epi32(_mm_max_epi32(intpart,gmin),gmax);
  _mm_storeu_si128 (&intpart, ixiyiz);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];
#else
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);
  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);
  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  __m128  uxuyuz= _mm_set_ps (ux, uy, uz, 0.0);
  __m128i intpart=_mm_set_epi32(ix,iy,iz,0);
#endif
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
  int Nm = (N+3)/4;
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
          mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[3*n+0] = ((float*)mgrad)[nd4*12 + 4*0 + nm4] * dxInv;
    grads[3*n+1] = ((float*)mgrad)[nd4*12 + 4*1 + nm4] * dyInv;
    grads[3*n+2] = ((float*)mgrad)[nd4*12 + 4*2 + nm4] * dzInv;
  }
}



void
eval_multi_UBspline_3d_s_vgl (const multi_UBspline_3d_s *spline,
                              float x, float y, float z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict lapl)
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
#if defined(HAVE_SSE41)
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
  //__m128i ixiyiz;
  //_mm_storeu_si128 (&ixiyiz, intpart);
  __m128i gmin = _mm_set_epi32(0,0,0,0);
  __m128i gmax = _mm_set_epi32(spline->x_grid.num-1, spline->y_grid.num-1,spline->z_grid.num-1,0);
  __m128i ixiyiz=_mm_min_epi32(_mm_max_epi32(intpart,gmin),gmax);
  _mm_storeu_si128 (&intpart, ixiyiz);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];
#else
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);
  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);
  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  __m128  uxuyuz= _mm_set_ps (ux, uy, uz, 0.0);
  __m128i intpart=_mm_set_epi32(ix,iy,iz,0);
#endif
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4,
         cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
         tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
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
  int Nm = (N+3)/4;
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
          mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[3*n+0] = ((float*)mgrad)[nd4*12 + 4*0 + nm4] * dxInv;
    grads[3*n+1] = ((float*)mgrad)[nd4*12 + 4*1 + nm4] * dyInv;
    grads[3*n+2] = ((float*)mgrad)[nd4*12 + 4*2 + nm4] * dzInv;
    lapl [n]     = (((float*)mlapl)[nd4*12 + 4*0 + nm4] * dxInv*dxInv +
                    ((float*)mlapl)[nd4*12 + 4*1 + nm4] * dyInv*dyInv +
                    ((float*)mlapl)[nd4*12 + 4*2 + nm4] * dzInv*dzInv);
  }
}



void
eval_multi_UBspline_3d_s_vgh (const multi_UBspline_3d_s *spline,
                              float x, float y, float z,
                              float* restrict vals,
                              float* restrict grads,
                              float* restrict hess)
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
#if defined(HAVE_SSE41)
  /// SSE mesh point determination
  __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
  __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start,
                                 spline->z_grid.start, 0.0);
  __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv,
                                 spline->z_grid.delta_inv, 0.0);
  xyz = _mm_sub_ps (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
  __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
  //__m128i ixiyiz;
  //_mm_storeu_si128 (&ixiyiz, intpart);
  __m128i gmin = _mm_set_epi32(0,0,0,0);
  __m128i gmax = _mm_set_epi32(spline->x_grid.num-1, spline->y_grid.num-1,spline->z_grid.num-1,0);
  __m128i ixiyiz=_mm_min_epi32(_mm_max_epi32(intpart,gmin),gmax);
  _mm_storeu_si128 (&intpart, ixiyiz);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];
#else
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);
  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);
  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  __m128  uxuyuz= _mm_set_ps (ux, uy, uz, 0.0);
  __m128i intpart=_mm_set_epi32(ix,iy,iz,0);
#endif
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
  int Nm = (N+3)/4;
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
          mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
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
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[3*n+0] = ((float*)mgrad)[nd4*12 + 4*0 + nm4] * dxInv;
    grads[3*n+1] = ((float*)mgrad)[nd4*12 + 4*1 + nm4] * dyInv;
    grads[3*n+2] = ((float*)mgrad)[nd4*12 + 4*2 + nm4] * dzInv;
    hess [9*n+0]               = ((float*)mhess)[nd4*24 + 4*0 + nm4] * dxInv*dxInv;
    hess [9*n+1] = hess[9*n+3] = ((float*)mhess)[nd4*24 + 4*1 + nm4] * dxInv*dyInv;
    hess [9*n+2] = hess[9*n+6] = ((float*)mhess)[nd4*24 + 4*2 + nm4] * dxInv*dzInv;
    hess [9*n+4]               = ((float*)mhess)[nd4*24 + 4*3 + nm4] * dyInv*dyInv;
    hess [9*n+5] = hess[9*n+7] = ((float*)mhess)[nd4*24 + 4*4 + nm4] * dyInv*dzInv;
    hess [9*n+8]               = ((float*)mhess)[nd4*24 + 4*5 + nm4] * dzInv*dzInv;
  }
//   for (int n=0; n<N/4; n++) {
//     _mm_storeu_ps ((float*)&vals[n*4], mvals[n]);
//     for (int i=0; i<3; i++) {
//       _mm_store_ss  ((float*)&grads[12*n+i+0], mgrad[3*n+i]);
//       _mm_store_ss  ((float*)&grads[12*n+i+3],_mm_shuffle_ps (mgrad[3*n+i], mgrad[3*n+i],_MM_SHUFFLE(1,1,1,1)));
//       _mm_store_ss  ((float*)&grads[12*n+i+6],_mm_shuffle_ps (mgrad[3*n+i], mgrad[3*n+i],_MM_SHUFFLE(2,2,2,2)));
//       _mm_store_ss  ((float*)&grads[12*n+i+9],_mm_shuffle_ps (mgrad[3*n+i], mgrad[3*n+i],_MM_SHUFFLE(3,3,3,3)));
//     }
//     _mm_store_ss  ((float*)&hess[36*n+0+ 0],                mhess[6*n+0]);
//     _mm_store_ss  ((float*)&hess[36*n+0+ 9],_mm_shuffle_ps (mhess[6*n+0], mhess[6*n+0],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+0+18],_mm_shuffle_ps (mhess[6*n+0], mhess[6*n+0],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+0+27],_mm_shuffle_ps (mhess[6*n+0], mhess[6*n+0],_MM_SHUFFLE(3,3,3,3)));
//     _mm_store_ss  ((float*)&hess[36*n+1+ 0],                mhess[6*n+1]);
//     _mm_store_ss  ((float*)&hess[36*n+1+ 9],_mm_shuffle_ps (mhess[6*n+1], mhess[6*n+1],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+1+18],_mm_shuffle_ps (mhess[6*n+1], mhess[6*n+1],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+1+27],_mm_shuffle_ps (mhess[6*n+1], mhess[6*n+1],_MM_SHUFFLE(3,3,3,3)));
//     _mm_store_ss  ((float*)&hess[36*n+2+ 0],                mhess[6*n+2]);
//     _mm_store_ss  ((float*)&hess[36*n+2+ 9],_mm_shuffle_ps (mhess[6*n+2], mhess[6*n+2],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+2+18],_mm_shuffle_ps (mhess[6*n+2], mhess[6*n+2],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+2+27],_mm_shuffle_ps (mhess[6*n+2], mhess[6*n+2],_MM_SHUFFLE(3,3,3,3)));
//     _mm_store_ss  ((float*)&hess[36*n+4+ 0],                mhess[6*n+3]);
//     _mm_store_ss  ((float*)&hess[36*n+4+ 9],_mm_shuffle_ps (mhess[6*n+3], mhess[6*n+3],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+4+18],_mm_shuffle_ps (mhess[6*n+3], mhess[6*n+3],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+4+27],_mm_shuffle_ps (mhess[6*n+3], mhess[6*n+3],_MM_SHUFFLE(3,3,3,3)));
//     _mm_store_ss  ((float*)&hess[36*n+5+ 0],                mhess[6*n+4]);
//     _mm_store_ss  ((float*)&hess[36*n+5+ 9],_mm_shuffle_ps (mhess[6*n+4], mhess[6*n+4],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+5+18],_mm_shuffle_ps (mhess[6*n+4], mhess[6*n+4],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+5+27],_mm_shuffle_ps (mhess[6*n+4], mhess[6*n+4],_MM_SHUFFLE(3,3,3,3)));
//     _mm_store_ss  ((float*)&hess[36*n+8+ 0],                mhess[6*n+5]);
//     _mm_store_ss  ((float*)&hess[36*n+8+ 9],_mm_shuffle_ps (mhess[6*n+5], mhess[6*n+5],_MM_SHUFFLE(1,1,1,1)));
//     _mm_store_ss  ((float*)&hess[36*n+8+18],_mm_shuffle_ps (mhess[6*n+5], mhess[6*n+5],_MM_SHUFFLE(2,2,2,2)));
//     _mm_store_ss  ((float*)&hess[36*n+8+27],_mm_shuffle_ps (mhess[6*n+5], mhess[6*n+5],_MM_SHUFFLE(3,3,3,3)));
//   }
//   // Store remainders
//   int mlast = Nm-1;
//   int nlast = 4*mlast;
//   for (int n=0; n<(N%4); n++) {
//     vals[nlast+n] = ((float*)mvals)[nlast+n];
//     for (int i=0; i<3; i++)
//       grads[3*(nlast+n)+i] = ((float*)mgrad)[3*nlast+4*i+n];
//     hess[9*(nlast+n)+0] = ((float*)mhess)[6*nlast+4*0+n];
//     hess[9*(nlast+n)+1] = ((float*)mhess)[6*nlast+4*1+n];
//     hess[9*(nlast+n)+2] = ((float*)mhess)[6*nlast+4*2+n];
//     hess[9*(nlast+n)+4] = ((float*)mhess)[6*nlast+4*3+n];
//     hess[9*(nlast+n)+5] = ((float*)mhess)[6*nlast+4*4+n];
//     hess[9*(nlast+n)+8] = ((float*)mhess)[6*nlast+4*5+n];
//   }
//   for (int n=0; n<N; n++) {
//     grads[3*n+0] *= dxInv;
//     grads[3*n+1] *= dyInv;
//     grads[3*n+2] *= dzInv;
//     hess [9*n+0] *= dxInv * dxInv;
//     hess [9*n+1] *= dxInv * dyInv;
//     hess [9*n+2] *= dxInv * dzInv;
//     hess [9*n+4] *= dyInv * dyInv;
//     hess [9*n+5] *= dyInv * dzInv;
//     hess [9*n+8] *= dzInv * dzInv;
//     hess [9*n+3] = hess[9*n+1];
//     hess [9*n+6] = hess[9*n+2];
//     hess [9*n+7] = hess[9*n+5];
//   }
}

inline void
eval_multi_UBspline_3d_s_vghgh (const multi_UBspline_3d_s *spline,
                                float x, float y, float z,
                                float* restrict vals,
                                float* restrict grads,
                                float* restrict hess,
                                float* restrict gradhess)
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
  _mm_prefetch ((const char*)  &A_s[12],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[14],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_s[15],_MM_HINT_T0);
#if defined(HAVE_SSE41)
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
  //__m128i ixiyiz;
  //_mm_storeu_si128 (&ixiyiz, intpart);
  __m128i gmin = _mm_set_epi32(0,0,0,0);
  __m128i gmax = _mm_set_epi32(spline->x_grid.num-1, spline->y_grid.num-1,spline->z_grid.num-1,0);
  __m128i ixiyiz=_mm_min_epi32(_mm_max_epi32(intpart,gmin),gmax);
  _mm_storeu_si128 (&intpart, ixiyiz);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix = ((int *)&ixiyiz)[3];
  int iy = ((int *)&ixiyiz)[2];
  int iz = ((int *)&ixiyiz)[1];
#else
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);
  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);
  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  __m128  uxuyuz= _mm_set_ps (ux, uy, uz, 0.0);
  __m128i intpart=_mm_set_epi32(ix,iy,iz,0);
#endif
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
  __m128 a4, b4, c4, da4, db4, dc4, d2a4, d2b4, d2c4, d3a4, d3b4, d3c4;
  // x-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpx, d2a4);
  _MM_MATVEC4_PS (A_s[12], A_s[13], A_s[14], A_s[15], tpx, d3a4);
  // y-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpy, d2b4);
  _MM_MATVEC4_PS (A_s[12], A_s[13], A_s[14], A_s[15], tpy, d3b4);
  // z-dependent vectors
  _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
  _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc4);
  _MM_MATVEC4_PS (A_s[ 8], A_s[ 9], A_s[10], A_s[11], tpz, d2c4);
  _MM_MATVEC4_PS (A_s[12], A_s[13], A_s[14], A_s[15], tpz, d3c4);
  __m128 a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4], d3a[4], d3b[4], d3c[4];
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
  tmp=_mm_unpacklo_ps(d3a4, d3a4);
  d3a[0]=_mm_unpacklo_ps(tmp, tmp);
  d3a[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d3a4, d3a4);
  d3a[2]=_mm_unpacklo_ps(tmp, tmp);
  d3a[3]=_mm_unpackhi_ps(tmp, tmp);
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
  tmp=_mm_unpacklo_ps(d3b4, d3b4);
  d3b[0]=_mm_unpacklo_ps(tmp, tmp);
  d3b[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d3b4, d3b4);
  d3b[2]=_mm_unpacklo_ps(tmp, tmp);
  d3b[3]=_mm_unpackhi_ps(tmp, tmp);
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
  tmp=_mm_unpacklo_ps(d3c4, d3c4);
  d3c[0]=_mm_unpacklo_ps(tmp, tmp);
  d3c[1]=_mm_unpackhi_ps(tmp, tmp);
  tmp=_mm_unpackhi_ps(d3c4, d3c4);
  d3c[2]=_mm_unpacklo_ps(tmp, tmp);
  d3c[3]=_mm_unpackhi_ps(tmp, tmp);
  int N = spline->num_splines;
  int Nm = (N+3)/4;
  __m128 mvals[Nm], mgrad[3*Nm], mhess[6*Nm], mgradhess[10*Nm];
  // Zero out values;
  __m128 mzero = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
  for (int n=0; n<Nm; n++)
    mvals[n] = _mm_setzero_ps();
  for (int n=0; n<3*Nm; n++)
    mgrad[n] = _mm_setzero_ps();
  for (int n=0; n<6*Nm; n++)
    mhess[n] = _mm_setzero_ps();
  for (int n=0; n<10*Nm; n++)
    mgradhess[n] = _mm_setzero_ps();
  // Main compute loop
  __m128 abc, dabc[3], d2abc[6],  d3abc[10];
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
        d3abc[0] = _mm_mul_ps (d3a[i], _mm_mul_ps(  b[j],  c[k]));
        d3abc[1] = _mm_mul_ps (d2a[i], _mm_mul_ps( db[j],  c[k]));
        d3abc[2] = _mm_mul_ps (d2a[i], _mm_mul_ps(  b[j],  dc[k]));
        d3abc[3] = _mm_mul_ps ( da[i], _mm_mul_ps(d2b[j],  c[k]));
        d3abc[4] = _mm_mul_ps ( da[i], _mm_mul_ps( db[j], dc[k]));
        d3abc[5] = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],d2c[k]));
        d3abc[6] = _mm_mul_ps (  a[i], _mm_mul_ps(d3b[j],  c[k]));
        d3abc[7] = _mm_mul_ps (  a[i], _mm_mul_ps(d2b[j], dc[k]));
        d3abc[8] = _mm_mul_ps (  a[i], _mm_mul_ps( db[j],d2c[k]));
        d3abc[9] = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],d3c[k]));
        __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nm; n++)
        {
          mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
          mgrad[3*n+0] = _mm_add_ps (mgrad[3*n+0], _mm_mul_ps( dabc[0], coefs[n]));
          mgrad[3*n+1] = _mm_add_ps (mgrad[3*n+1], _mm_mul_ps( dabc[1], coefs[n]));
          mgrad[3*n+2] = _mm_add_ps (mgrad[3*n+2], _mm_mul_ps( dabc[2], coefs[n]));
          mhess[6*n+0] = _mm_add_ps (mhess[6*n+0], _mm_mul_ps(d2abc[0], coefs[n]));
          mhess[6*n+1] = _mm_add_ps (mhess[6*n+1], _mm_mul_ps(d2abc[1], coefs[n]));
          mhess[6*n+2] = _mm_add_ps (mhess[6*n+2], _mm_mul_ps(d2abc[2], coefs[n]));
          mhess[6*n+3] = _mm_add_ps (mhess[6*n+3], _mm_mul_ps(d2abc[3], coefs[n]));
          mhess[6*n+4] = _mm_add_ps (mhess[6*n+4], _mm_mul_ps(d2abc[4], coefs[n]));
          mhess[6*n+5] = _mm_add_ps (mhess[6*n+5], _mm_mul_ps(d2abc[5], coefs[n]));
          mgradhess[10*n+0] = _mm_add_ps (mgradhess[10*n+0], _mm_mul_ps(d3abc[0], coefs[n]));
          mgradhess[10*n+1] = _mm_add_ps (mgradhess[10*n+1], _mm_mul_ps(d3abc[1], coefs[n]));
          mgradhess[10*n+2] = _mm_add_ps (mgradhess[10*n+2], _mm_mul_ps(d3abc[2], coefs[n]));
          mgradhess[10*n+3] = _mm_add_ps (mgradhess[10*n+3], _mm_mul_ps(d3abc[3], coefs[n]));
          mgradhess[10*n+4] = _mm_add_ps (mgradhess[10*n+4], _mm_mul_ps(d3abc[4], coefs[n]));
          mgradhess[10*n+5] = _mm_add_ps (mgradhess[10*n+5], _mm_mul_ps(d3abc[5], coefs[n]));
          mgradhess[10*n+6] = _mm_add_ps (mgradhess[10*n+6], _mm_mul_ps(d3abc[6], coefs[n]));
          mgradhess[10*n+7] = _mm_add_ps (mgradhess[10*n+7], _mm_mul_ps(d3abc[7], coefs[n]));
          mgradhess[10*n+8] = _mm_add_ps (mgradhess[10*n+8], _mm_mul_ps(d3abc[8], coefs[n]));
          mgradhess[10*n+9] = _mm_add_ps (mgradhess[10*n+9], _mm_mul_ps(d3abc[9], coefs[n]));
        }
      }
  // Now, store results back
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N; n++)
  {
    int nd4 = n>>2;
    int nm4 = n & 3;
    vals[n] = ((float*)mvals)[n];
    grads[3*n+0] = ((float*)mgrad)[nd4*12 + 4*0 + nm4] * dxInv;
    grads[3*n+1] = ((float*)mgrad)[nd4*12 + 4*1 + nm4] * dyInv;
    grads[3*n+2] = ((float*)mgrad)[nd4*12 + 4*2 + nm4] * dzInv;
    hess [9*n+0]               = ((float*)mhess)[nd4*24 + 4*0 + nm4] * dxInv*dxInv;
    hess [9*n+1] = hess[9*n+3] = ((float*)mhess)[nd4*24 + 4*1 + nm4] * dxInv*dyInv;
    hess [9*n+2] = hess[9*n+6] = ((float*)mhess)[nd4*24 + 4*2 + nm4] * dxInv*dzInv;
    hess [9*n+4]               = ((float*)mhess)[nd4*24 + 4*3 + nm4] * dyInv*dyInv;
    hess [9*n+5] = hess[9*n+7] = ((float*)mhess)[nd4*24 + 4*4 + nm4] * dyInv*dzInv;
    hess [9*n+8]               = ((float*)mhess)[nd4*24 + 4*5 + nm4] * dzInv*dzInv;
    gradhess [27*n+0 ] =   ((float*)mgradhess)[nd4*40 + 4*0 + nm4] * dxInv*dxInv*dxInv;
    gradhess [27*n+1 ] =   ((float*)mgradhess)[nd4*40 + 4*1 + nm4] * dxInv*dxInv*dyInv;
    gradhess [27*n+2 ] =   ((float*)mgradhess)[nd4*40 + 4*2 + nm4] * dxInv*dxInv*dzInv;
    gradhess [27*n+4 ] =   ((float*)mgradhess)[nd4*40 + 4*3 + nm4] * dxInv*dyInv*dyInv;
    gradhess [27*n+5 ] =   ((float*)mgradhess)[nd4*40 + 4*4 + nm4] * dxInv*dyInv*dzInv;
    gradhess [27*n+8 ] =   ((float*)mgradhess)[nd4*40 + 4*5 + nm4] * dxInv*dzInv*dzInv;
    gradhess [27*n+13] =   ((float*)mgradhess)[nd4*40 + 4*6 + nm4] * dyInv*dyInv*dyInv;
    gradhess [27*n+14] =   ((float*)mgradhess)[nd4*40 + 4*7 + nm4] * dyInv*dyInv*dzInv;
    gradhess [27*n+17] =   ((float*)mgradhess)[nd4*40 + 4*8 + nm4] * dyInv*dzInv*dzInv;
    gradhess [27*n+26] =   ((float*)mgradhess)[nd4*40 + 4*9 + nm4] * dzInv*dzInv*dzInv;
    // Copy gradhess elements into rest of tensor
    gradhess [27*n+9  ] = gradhess [27*n+3  ] = gradhess [27*n+1 ];
    gradhess [27*n+18 ] = gradhess [27*n+6  ] = gradhess [27*n+2 ];
    gradhess [27*n+22 ] = gradhess [27*n+16 ] = gradhess [27*n+14];
    gradhess [27*n+12 ] = gradhess [27*n+10 ] = gradhess [27*n+4 ];
    gradhess [27*n+24 ] = gradhess [27*n+20 ] = gradhess [27*n+8 ];
    gradhess [27*n+25 ] = gradhess [27*n+23 ] = gradhess [27*n+17];
    gradhess [27*n+21 ] = gradhess [27*n+19 ] = gradhess [27*n+15] = gradhess [27*n+11 ] = gradhess [27*n+7 ] = gradhess [27*n+5];
  }
}
