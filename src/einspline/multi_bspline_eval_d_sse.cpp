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

#include "config.h"

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#include <math.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"
#include "multi_bspline_eval_d.h"

extern __m128d *restrict A_d;
extern double *restrict Ad, *restrict dAd, *restrict d2Ad, *restrict d3Ad;

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

/*********************************************************/
/* 1D double-precision, real evaulation functions        */
/*********************************************************/
#include <stdio.h>
void
eval_multi_UBspline_1d_d (const multi_UBspline_1d_d *spline,
                          double x,
                          double* restrict vals)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  double tpx[4], a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  double* restrict coefs = spline->coefs;
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  // fprintf (stderr, "ux = %12.8f\n", ux);
  // fprintf (stderr, "ipart = %ix  tx = %12.7f\n", ix, tx);
  // fprintf (stderr, "a[0] = %1.8e\n", a[0]);
  // fprintf (stderr, "a[1] = %1.8e\n", a[1]);
  // fprintf (stderr, "a[2] = %1.8e\n", a[2]);
  // fprintf (stderr, "a[3] = %1.8e\n", a[3]);
  // fprintf (stderr, "tpx[0] = %1.8e\n", tpx[0]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
    vals[n]  = 0.0;
  for (int i=0; i<4; i++)
  {
    double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
      vals[n]  +=   a[i] * coefs[n];
  }
}



void
eval_multi_UBspline_1d_d_vg (const multi_UBspline_1d_d *spline,
                             double x,
                             double* restrict vals,
                             double* restrict grads)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  double tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  double* restrict coefs = spline->coefs;
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
    {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
    }
  }
  double dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
    grads[n] *= dxInv;
}


void
eval_multi_UBspline_1d_d_vgl (const multi_UBspline_1d_d *spline,
                              double x,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  double tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  double* restrict coefs = spline->coefs;
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
    {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
      lapl[n]  += d2a[i] * coefs[n];
    }
  }
  double dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
    grads[n] *= dxInv;
    lapl [n] *= dxInv*dxInv;
  }
}



void
eval_multi_UBspline_1d_d_vgh (const multi_UBspline_1d_d *spline,
                              double x,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess)
{
  eval_multi_UBspline_1d_d_vgl (spline, x, vals, grads, hess);
}

/*********************************************************/
/* 2D double-precision, real evaulation functions        */
/*********************************************************/
void
eval_multi_UBspline_2d_d(const multi_UBspline_2d_d *spline,
                         double x, double y,
                         double* restrict vals)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01  ,   b01,   a23,    b23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh];
  for (int n=0; n<Nh; n++)
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
  __m128d a[4], b[4], da[4], db[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      __m128d ab;
      ab         = _mm_mul_pd(a[i], b[j]);
      __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nh; n++)
        mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
    }
  for (int n=0; n<N/2; n++)
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
  if (N&1)
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
}


void
eval_multi_UBspline_2d_d_vg (const multi_UBspline_2d_d *spline,
                             double x, double y,
                             double* restrict vals,
                             double* restrict grads)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01  ,   b01,   a23,    b23,
          da01 ,  db01,  da23,   db23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[2*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<2; i++)
      mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i], mgrads[2*n+i]);
  }
  __m128d a[4], b[4], da[4], db[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      __m128d ab, d_ab[2];
      ab         = _mm_mul_pd(a[i], b[j]);
      d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
      d_ab[1]    = _mm_mul_pd( a[i],  db[j]);
      __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nh; n++)
      {
        mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
        mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
        mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+4*n+0), mgrads[2*n+0]);
    _mm_storeh_pd((double*)(grads+4*n+2), mgrads[2*n+0]);
    _mm_storel_pd((double*)(grads+4*n+1), mgrads[2*n+1]);
    _mm_storeh_pd((double*)(grads+4*n+3), mgrads[2*n+1]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+2*(N-1)+0),mgrads[2*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+2*(N-1)+1),mgrads[2*(Nh-1)+1]);
  }
  for (int n=0; n<N; n++)
  {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
  }
}



void
eval_multi_UBspline_2d_d_vgl (const multi_UBspline_2d_d *spline,
                              double x, double y,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01  ,   b01,   a23,    b23,
          da01 ,  db01,  da23,   db23,
          d2a01, d2b01, d2a23,  d2b23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[2*Nh], mlapl[2*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<2; i++)
    {
      mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i], mgrads[2*n+i]);
      mlapl [2*n+i] = _mm_sub_pd (mlapl [2*n+i], mlapl [2*n+i]);
    }
  }
  __m128d a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      __m128d ab, d_ab[2], d2_ab[2];
      ab         = _mm_mul_pd(a[i], b[j]);
      d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
      d_ab[1]    = _mm_mul_pd( a[i],  db[j]);
      d2_ab[0]   = _mm_mul_pd(d2a[i],   b[j]);
      d2_ab[1]   = _mm_mul_pd(  a[i], d2b[j]);
      __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nh; n++)
      {
        mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
        mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
        mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
        mlapl[2*n+0]  = _mm_add_pd (mlapl[2*n+0],  _mm_mul_pd (d2_ab[0], coefs[n]));
        mlapl[2*n+1]  = _mm_add_pd (mlapl[2*n+1],  _mm_mul_pd (d2_ab[1], coefs[n]));
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double lapl2[2*N];
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+4*n+0), mgrads[2*n+0]);
    _mm_storeh_pd((double*)(grads+4*n+2), mgrads[2*n+0]);
    _mm_storel_pd((double*)(grads+4*n+1), mgrads[2*n+1]);
    _mm_storeh_pd((double*)(grads+4*n+3), mgrads[2*n+1]);
    _mm_storel_pd((double*)(lapl2+4*n+0), mlapl [2*n+0]);
    _mm_storeh_pd((double*)(lapl2+4*n+2), mlapl [2*n+0]);
    _mm_storel_pd((double*)(lapl2+4*n+1), mlapl [2*n+1]);
    _mm_storeh_pd((double*)(lapl2+4*n+3), mlapl [2*n+1]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+2*(N-1)+0),mgrads[2*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+2*(N-1)+1),mgrads[2*(Nh-1)+1]);
    _mm_storel_pd((double*)(lapl2+2*(N-1)+0),  mlapl [2*(Nh-1)+0]);
    _mm_storel_pd((double*)(lapl2+2*(N-1)+1),  mlapl [2*(Nh-1)+1]);
  }
  for (int n=0; n<N; n++)
  {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    lapl2[2*n+0] *= dxInv*dxInv;
    lapl2[2*n+1] *= dyInv*dyInv;
    lapl [n] = lapl2[2*n+0] + lapl2[2*n+1];
  }
}



void
eval_multi_UBspline_2d_d_vgh (const multi_UBspline_2d_d *spline,
                              double x, double y,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01  ,   b01,   a23,    b23,
          da01 ,  db01,  da23,   db23,
          d2a01, d2b01, d2a23,  d2b23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[2*Nh], mhess[3*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<2; i++)
      mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i],mgrads[2*n+i]);
    for (int i=0; i<3; i++)
      mhess[3*n+i]  = _mm_sub_pd (mhess[3*n+i], mhess[3*n+i]);
  }
  __m128d a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      __m128d ab, d_ab[2], d2_ab[3];
      ab         = _mm_mul_pd(a[i], b[j]);
      d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
      d_ab[1]    = _mm_mul_pd( a[i],  db[j]);
      d2_ab[0]   = _mm_mul_pd(d2a[i],   b[j]);
      d2_ab[1]   = _mm_mul_pd( da[i],  db[j]);
      d2_ab[2]   = _mm_mul_pd(  a[i], d2b[j]);
      __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<Nh; n++)
      {
        mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
        mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
        mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
        mhess[3*n+0]  = _mm_add_pd (mhess[3*n+0],  _mm_mul_pd (d2_ab[0], coefs[n]));
        mhess[3*n+1]  = _mm_add_pd (mhess[3*n+1],  _mm_mul_pd (d2_ab[1], coefs[n]));
        mhess[3*n+2]  = _mm_add_pd (mhess[3*n+2],  _mm_mul_pd (d2_ab[2], coefs[n]));
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+4*n+0),mgrads[2*n+0]);
    _mm_storeh_pd((double*)(grads+4*n+2),mgrads[2*n+0]);
    _mm_storel_pd((double*)(grads+4*n+1),mgrads[2*n+1]);
    _mm_storeh_pd((double*)(grads+4*n+3),mgrads[2*n+1]);
    _mm_storel_pd((double*)(hess+8*n+0),  mhess [3*n+0]);
    _mm_storeh_pd((double*)(hess+8*n+4),  mhess [3*n+0]);
    _mm_storel_pd((double*)(hess+8*n+1),  mhess [3*n+1]);
    _mm_storeh_pd((double*)(hess+8*n+5),  mhess [3*n+1]);
    _mm_storel_pd((double*)(hess+8*n+3),  mhess [3*n+2]);
    _mm_storeh_pd((double*)(hess+8*n+7),  mhess [3*n+2]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+2*(N-1)+0),mgrads[2*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+2*(N-1)+1),mgrads[2*(Nh-1)+1]);
    _mm_storel_pd((double*)(hess+4*(N-1)+0),  mhess [3*(Nh-1)+0]);
    _mm_storel_pd((double*)(hess+4*(N-1)+1),  mhess [3*(Nh-1)+1]);
    _mm_storel_pd((double*)(hess+4*(N-1)+3),  mhess [3*(Nh-1)+2]);
  }
  for (int n=0; n<N; n++)
  {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    hess[4*n+0]  *= dxInv*dxInv;
    hess[4*n+1]  *= dxInv*dyInv;
    hess[4*n+3]  *= dyInv*dyInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess[4*n+2] = hess[4*n+1];
  }
}


/*********************************************************/
/* 3D double-precision, real evaulation functions        */
/*********************************************************/
void
eval_multi_UBspline_3d_d (const multi_UBspline_3d_d *spline,
                          double x, double y, double z,
                          double* restrict vals)
{
  _mm_prefetch ((const char*) &A_d[0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[7],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01, b01, c01, a23, b23, c23,
          tmp0, tmp1, r0, r1, i0, i1, val_r, val_i;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpx01, tpx23, tpx01, tpx23,   a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpy01, tpy23, tpy01, tpy23,   b23);
  // z-dependent vectors
  _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpz01, tpz23, tpz01, tpz23,   c01);
  _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpz01, tpz23, tpz01, tpz23,   c23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh];
  for (int n=0; n<Nh; n++)
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
  __m128d a[4], b[4], c[4];
  a[0]   = _mm_unpacklo_pd(a01,a01);
  a[1]   = _mm_unpackhi_pd(a01,a01);
  a[2]   = _mm_unpacklo_pd(a23,a23);
  a[3]   = _mm_unpackhi_pd(a23,a23);
  b[0]   = _mm_unpacklo_pd(b01,b01);
  b[1]   = _mm_unpackhi_pd(b01,b01);
  b[2]   = _mm_unpacklo_pd(b23,b23);
  b[3]   = _mm_unpackhi_pd(b23,b23);
  c[0]   = _mm_unpacklo_pd(c01,c01);
  c[1]   = _mm_unpackhi_pd(c01,c01);
  c[2]   = _mm_unpacklo_pd(c23,c23);
  c[3]   = _mm_unpackhi_pd(c23,c23);
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        __m128d abc = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
        __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nh; n++)
          mvals[n] = _mm_add_pd (mvals[n], _mm_mul_pd (abc, coefs[n]));
      }
  for (int n=0; n<N/2; n++)
    _mm_storeu_pd((vals+2*n),mvals[n]);
  if (N & 1)
    _mm_storel_pd(vals+N-1,mvals[N/2]);
}



void
eval_multi_UBspline_3d_d_vg (const multi_UBspline_3d_d *spline,
                             double x, double y, double z,
                             double* restrict vals,
                             double* restrict grads)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01  ,   b01,   c01,   a23,    b23,   c23,
          da01 ,  db01,  dc01,  da23,   db23,  dc23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  // z-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[3*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<3; i++)
      mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
  }
  __m128d a[4], b[4], c[4], da[4], db[4], dc[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  c[0]=_mm_unpacklo_pd(c01,c01);
  dc[0]=_mm_unpacklo_pd(dc01,dc01);
  c[1]=_mm_unpackhi_pd(c01,c01);
  dc[1]=_mm_unpackhi_pd(dc01,dc01);
  c[2]=_mm_unpacklo_pd(c23,c23);
  dc[2]=_mm_unpacklo_pd(dc23,dc23);
  c[3]=_mm_unpackhi_pd(c23,c23);
  dc[3]=_mm_unpackhi_pd(dc23,dc23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        __m128d abc, d_abc[3];
        abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
        d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
        d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
        d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);
        __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nh; n++)
        {
          mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
          mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
          mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
          mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
        }
      }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+6*n+0),mgrads[3*n+0]);
    _mm_storeh_pd((double*)(grads+6*n+3),mgrads[3*n+0]);
    _mm_storel_pd((double*)(grads+6*n+1),mgrads[3*n+1]);
    _mm_storeh_pd((double*)(grads+6*n+4),mgrads[3*n+1]);
    _mm_storel_pd((double*)(grads+6*n+2),mgrads[3*n+2]);
    _mm_storeh_pd((double*)(grads+6*n+5),mgrads[3*n+2]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+0),mgrads[3*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+3*(N-1)+1),mgrads[3*(Nh-1)+1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+2),mgrads[3*(Nh-1)+2]);
  }
  for (int n=0; n<N; n++)
  {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
  }
}



void
eval_multi_UBspline_3d_d_vgl (const multi_UBspline_3d_d *spline,
                              double x, double y, double z,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict lapl)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01  ,   b01,   c01,   a23,    b23,   c23,
          da01 ,  db01,  dc01,  da23,   db23,  dc23,
          d2a01, d2b01, d2c01, d2a23,  d2b23, d2c23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);
  // z-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpz01, tpz23, tpz01, tpz23, d2c01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpz01, tpz23, tpz01, tpz23, d2c23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[3*Nh], mlapl[3*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<3; i++)
    {
      mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
      mlapl[3*n+i]  = _mm_sub_pd ( mlapl[3*n+i], mlapl[3*n+i]);
    }
  }
  __m128d a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
  c[0]=_mm_unpacklo_pd(c01,c01);
  dc[0]=_mm_unpacklo_pd(dc01,dc01);
  d2c[0]=_mm_unpacklo_pd(d2c01,d2c01);
  c[1]=_mm_unpackhi_pd(c01,c01);
  dc[1]=_mm_unpackhi_pd(dc01,dc01);
  d2c[1]=_mm_unpackhi_pd(d2c01,d2c01);
  c[2]=_mm_unpacklo_pd(c23,c23);
  dc[2]=_mm_unpacklo_pd(dc23,dc23);
  d2c[2]=_mm_unpacklo_pd(d2c23,d2c23);
  c[3]=_mm_unpackhi_pd(c23,c23);
  dc[3]=_mm_unpackhi_pd(dc23,dc23);
  d2c[3]=_mm_unpackhi_pd(d2c23,d2c23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        __m128d abc, d_abc[3], d2_abc[3];
        abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
        d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
        d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
        d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);
        d2_abc[0]   = _mm_mul_pd (_mm_mul_pd(d2a[i],   b[j]),   c[k]);
        d2_abc[1]   = _mm_mul_pd (_mm_mul_pd(  a[i], d2b[j]),   c[k]);
        d2_abc[2]   = _mm_mul_pd (_mm_mul_pd(  a[i],   b[j]), d2c[k]);
        __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nh; n++)
        {
          mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
          mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
          mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
          mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
          mlapl[3*n+0]  = _mm_add_pd (mlapl [3*n+0], _mm_mul_pd (d2_abc[0], coefs[n]));
          mlapl[3*n+1]  = _mm_add_pd (mlapl [3*n+1], _mm_mul_pd (d2_abc[1], coefs[n]));
          mlapl[3*n+2]  = _mm_add_pd (mlapl [3*n+2], _mm_mul_pd (d2_abc[2], coefs[n]));
        }
      }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  double lapl3[3*N];
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+6*n+0),mgrads[3*n+0]);
    _mm_storeh_pd((double*)(grads+6*n+3),mgrads[3*n+0]);
    _mm_storel_pd((double*)(grads+6*n+1),mgrads[3*n+1]);
    _mm_storeh_pd((double*)(grads+6*n+4),mgrads[3*n+1]);
    _mm_storel_pd((double*)(grads+6*n+2),mgrads[3*n+2]);
    _mm_storeh_pd((double*)(grads+6*n+5),mgrads[3*n+2]);
    _mm_storel_pd((double*)(lapl3+6*n+0), mlapl [3*n+0]);
    _mm_storeh_pd((double*)(lapl3+6*n+3), mlapl [3*n+0]);
    _mm_storel_pd((double*)(lapl3+6*n+1), mlapl [3*n+1]);
    _mm_storeh_pd((double*)(lapl3+6*n+4), mlapl [3*n+1]);
    _mm_storel_pd((double*)(lapl3+6*n+2), mlapl [3*n+2]);
    _mm_storeh_pd((double*)(lapl3+6*n+5), mlapl [3*n+2]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+0),mgrads[3*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+3*(N-1)+1),mgrads[3*(Nh-1)+1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+2),mgrads[3*(Nh-1)+2]);
    _mm_storel_pd((double*)(lapl3+3*(N-1)+0),  mlapl [3*(Nh-1)+0]);
    _mm_storel_pd((double*)(lapl3+3*(N-1)+1),  mlapl [3*(Nh-1)+1]);
    _mm_storel_pd((double*)(lapl3+3*(N-1)+2),  mlapl [3*(Nh-1)+2]);
  }
  for (int n=0; n<N; n++)
  {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    lapl3[3*n+0] *= dxInv*dxInv;
    lapl3[3*n+1] *= dyInv*dyInv;
    lapl3[3*n+2] *= dzInv*dzInv;
    lapl[n] = lapl3[3*n+0] + lapl3[3*n+1] + lapl3[3*n+2];
  }
}


void
eval_multi_UBspline_3d_d_vgh (const multi_UBspline_3d_d *spline,
                              double x, double y, double z,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01  ,   b01,   c01,   a23,    b23,   c23,
          da01 ,  db01,  dc01,  da23,   db23,  dc23,
          d2a01, d2b01, d2c01, d2a23,  d2b23, d2c23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);
  // z-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpz01, tpz23, tpz01, tpz23, d2c01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpz01, tpz23, tpz01, tpz23, d2c23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[3*Nh], mhess[6*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<3; i++)
      mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
    for (int i=0; i<6; i++)
      mhess[6*n+i]  = _mm_sub_pd (mhess[6*n+i], mhess[6*n+i]);
  }
  __m128d a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
  c[0]=_mm_unpacklo_pd(c01,c01);
  dc[0]=_mm_unpacklo_pd(dc01,dc01);
  d2c[0]=_mm_unpacklo_pd(d2c01,d2c01);
  c[1]=_mm_unpackhi_pd(c01,c01);
  dc[1]=_mm_unpackhi_pd(dc01,dc01);
  d2c[1]=_mm_unpackhi_pd(d2c01,d2c01);
  c[2]=_mm_unpacklo_pd(c23,c23);
  dc[2]=_mm_unpacklo_pd(dc23,dc23);
  d2c[2]=_mm_unpacklo_pd(d2c23,d2c23);
  c[3]=_mm_unpackhi_pd(c23,c23);
  dc[3]=_mm_unpackhi_pd(dc23,dc23);
  d2c[3]=_mm_unpackhi_pd(d2c23,d2c23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        __m128d abc, d_abc[3], d2_abc[6];
        abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
        d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
        d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
        d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);
        d2_abc[0]   = _mm_mul_pd (_mm_mul_pd(d2a[i],   b[j]),   c[k]);
        d2_abc[1]   = _mm_mul_pd (_mm_mul_pd( da[i],  db[j]),   c[k]);
        d2_abc[2]   = _mm_mul_pd (_mm_mul_pd( da[i],   b[j]),  dc[k]);
        d2_abc[3]   = _mm_mul_pd (_mm_mul_pd(  a[i], d2b[j]),   c[k]);
        d2_abc[4]   = _mm_mul_pd (_mm_mul_pd(  a[i],  db[j]),  dc[k]);
        d2_abc[5]   = _mm_mul_pd (_mm_mul_pd(  a[i],   b[j]), d2c[k]);
        __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nh; n++)
        {
          mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
          mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
          mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
          mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
          mhess[6*n+0]  = _mm_add_pd (mhess[6*n+0],  _mm_mul_pd (d2_abc[0], coefs[n]));
          mhess[6*n+1]  = _mm_add_pd (mhess[6*n+1],  _mm_mul_pd (d2_abc[1], coefs[n]));
          mhess[6*n+2]  = _mm_add_pd (mhess[6*n+2],  _mm_mul_pd (d2_abc[2], coefs[n]));
          mhess[6*n+3]  = _mm_add_pd (mhess[6*n+3],  _mm_mul_pd (d2_abc[3], coefs[n]));
          mhess[6*n+4]  = _mm_add_pd (mhess[6*n+4],  _mm_mul_pd (d2_abc[4], coefs[n]));
          mhess[6*n+5]  = _mm_add_pd (mhess[6*n+5],  _mm_mul_pd (d2_abc[5], coefs[n]));
        }
      }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+6*n+0),mgrads[3*n+0]);
    _mm_storeh_pd((double*)(grads+6*n+3),mgrads[3*n+0]);
    _mm_storel_pd((double*)(grads+6*n+1),mgrads[3*n+1]);
    _mm_storeh_pd((double*)(grads+6*n+4),mgrads[3*n+1]);
    _mm_storel_pd((double*)(grads+6*n+2),mgrads[3*n+2]);
    _mm_storeh_pd((double*)(grads+6*n+5),mgrads[3*n+2]);
    _mm_storel_pd((double*)(hess+18*n+0),  mhess [6*n+0]);
    _mm_storeh_pd((double*)(hess+18*n+9),  mhess [6*n+0]);
    _mm_storel_pd((double*)(hess+18*n+1),  mhess [6*n+1]);
    _mm_storeh_pd((double*)(hess+18*n+10), mhess [6*n+1]);
    _mm_storel_pd((double*)(hess+18*n+2),  mhess [6*n+2]);
    _mm_storeh_pd((double*)(hess+18*n+11), mhess [6*n+2]);
    _mm_storel_pd((double*)(hess+18*n+4),  mhess [6*n+3]);
    _mm_storeh_pd((double*)(hess+18*n+13), mhess [6*n+3]);
    _mm_storel_pd((double*)(hess+18*n+5),  mhess [6*n+4]);
    _mm_storeh_pd((double*)(hess+18*n+14), mhess [6*n+4]);
    _mm_storel_pd((double*)(hess+18*n+8),  mhess [6*n+5]);
    _mm_storeh_pd((double*)(hess+18*n+17), mhess [6*n+5]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+0),mgrads[3*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+3*(N-1)+1),mgrads[3*(Nh-1)+1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+2),mgrads[3*(Nh-1)+2]);
    _mm_storel_pd((double*)(hess+9*(N-1)+0),  mhess [6*(Nh-1)+0]);
    _mm_storel_pd((double*)(hess+9*(N-1)+1),  mhess [6*(Nh-1)+1]);
    _mm_storel_pd((double*)(hess+9*(N-1)+2),  mhess [6*(Nh-1)+2]);
    _mm_storel_pd((double*)(hess+9*(N-1)+4),  mhess [6*(Nh-1)+3]);
    _mm_storel_pd((double*)(hess+9*(N-1)+5),  mhess [6*(Nh-1)+4]);
    _mm_storel_pd((double*)(hess+9*(N-1)+8),  mhess [6*(Nh-1)+5]);
  }
  for (int n=0; n<N; n++)
  {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess[9*n+0]  *= dxInv*dxInv;
    hess[9*n+4]  *= dyInv*dyInv;
    hess[9*n+8]  *= dzInv*dzInv;
    hess[9*n+1]  *= dxInv*dyInv;
    hess[9*n+2]  *= dxInv*dzInv;
    hess[9*n+5]  *= dyInv*dzInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess[9*n+3] = hess[9*n+1];
    hess[9*n+6] = hess[9*n+2];
    hess[9*n+7] = hess[9*n+5];
  }
}


void
eval_multi_UBspline_3d_d_vghgh (const multi_UBspline_3d_d *spline,
                                double x, double y, double z,
                                double* restrict vals,
                                double* restrict grads,
                                double* restrict hess,
                                double* restrict gradhess)
{
  _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[24],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[25],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[26],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[27],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[28],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[29],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[30],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[31],_MM_HINT_T0);
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  int N  = spline->num_splines;
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01  ,   b01,   c01,   a23,    b23,   c23,
          da01 ,  db01,  dc01,  da23,   db23,  dc23,
          d2a01, d2b01, d2c01, d2a23,  d2b23, d2c23,
          d3a01, d3b01, d3c01, d3a23,  d3b23, d3c23;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
  tpz23 = _mm_set_pd (tz, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);
  _MM_DDOT4_PD (A_d[24], A_d[25], A_d[26], A_d[27], tpx01, tpx23, tpx01, tpx23, d3a01);
  _MM_DDOT4_PD (A_d[28], A_d[29], A_d[30], A_d[31], tpx01, tpx23, tpx01, tpx23, d3a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);
  _MM_DDOT4_PD (A_d[24], A_d[25], A_d[26], A_d[27], tpy01, tpy23, tpy01, tpy23, d3b01);
  _MM_DDOT4_PD (A_d[28], A_d[29], A_d[30], A_d[31], tpy01, tpy23, tpy01, tpy23, d3b23);
  // z-dependent vectors
  _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
  _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
  _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
  _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
  _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpz01, tpz23, tpz01, tpz23, d2c01);
  _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpz01, tpz23, tpz01, tpz23, d2c23);
  _MM_DDOT4_PD (A_d[24], A_d[25], A_d[26], A_d[27], tpz01, tpz23, tpz01, tpz23, d3c01);
  _MM_DDOT4_PD (A_d[28], A_d[29], A_d[30], A_d[31], tpz01, tpz23, tpz01, tpz23, d3c23);
  // Zero-out values
  int Nh = (N+1)/2;
  __m128d mvals[Nh], mgrads[3*Nh], mhess[6*Nh], mgradhess[10*Nh];
  for (int n=0; n<Nh; n++)
  {
    mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
    for (int i=0; i<3; i++)
      mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
    for (int i=0; i<6; i++)
      mhess[6*n+i]  = _mm_sub_pd (mhess[6*n+i], mhess[6*n+i]);
    for (int i=0; i<10; i++)
      mgradhess[10*n+i]  = _mm_sub_pd (mgradhess[10*n+i], mgradhess[10*n+i]);
  }
  __m128d a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4], d3a[4], d3b[4], d3c[4];
  a[0]=_mm_unpacklo_pd(a01,a01);
  da[0]=_mm_unpacklo_pd(da01,da01);
  d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
  d3a[0]=_mm_unpacklo_pd(d3a01,d3a01);
  a[1]=_mm_unpackhi_pd(a01,a01);
  da[1]=_mm_unpackhi_pd(da01,da01);
  d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
  d3a[1]=_mm_unpackhi_pd(d3a01,d3a01);
  a[2]=_mm_unpacklo_pd(a23,a23);
  da[2]=_mm_unpacklo_pd(da23,da23);
  d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
  d3a[2]=_mm_unpacklo_pd(d3a23,d3a23);
  a[3]=_mm_unpackhi_pd(a23,a23);
  da[3]=_mm_unpackhi_pd(da23,da23);
  d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
  d3a[3]=_mm_unpackhi_pd(d3a23,d3a23);
  b[0]=_mm_unpacklo_pd(b01,b01);
  db[0]=_mm_unpacklo_pd(db01,db01);
  d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
  d3b[0]=_mm_unpacklo_pd(d3b01,d3b01);
  b[1]=_mm_unpackhi_pd(b01,b01);
  db[1]=_mm_unpackhi_pd(db01,db01);
  d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
  d3b[1]=_mm_unpackhi_pd(d3b01,d3b01);
  b[2]=_mm_unpacklo_pd(b23,b23);
  db[2]=_mm_unpacklo_pd(db23,db23);
  d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
  d3b[2]=_mm_unpacklo_pd(d3b23,d3b23);
  b[3]=_mm_unpackhi_pd(b23,b23);
  db[3]=_mm_unpackhi_pd(db23,db23);
  d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
  d3b[3]=_mm_unpackhi_pd(d3b23,d3b23);
  c[0]=_mm_unpacklo_pd(c01,c01);
  dc[0]=_mm_unpacklo_pd(dc01,dc01);
  d2c[0]=_mm_unpacklo_pd(d2c01,d2c01);
  d3c[0]=_mm_unpacklo_pd(d3c01,d3c01);
  c[1]=_mm_unpackhi_pd(c01,c01);
  dc[1]=_mm_unpackhi_pd(dc01,dc01);
  d2c[1]=_mm_unpackhi_pd(d2c01,d2c01);
  d3c[1]=_mm_unpackhi_pd(d3c01,d3c01);
  c[2]=_mm_unpacklo_pd(c23,c23);
  dc[2]=_mm_unpacklo_pd(dc23,dc23);
  d2c[2]=_mm_unpacklo_pd(d2c23,d2c23);
  d3c[2]=_mm_unpacklo_pd(d3c23,d3c23);
  c[3]=_mm_unpackhi_pd(c23,c23);
  dc[3]=_mm_unpackhi_pd(dc23,dc23);
  d2c[3]=_mm_unpackhi_pd(d2c23,d2c23);
  d3c[3]=_mm_unpackhi_pd(d3c23,d3c23);
  // Main computation loop
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        __m128d abc, d_abc[3], d2_abc[6], d3_abc[10];
        abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
        d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
        d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
        d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);
        d2_abc[0]   = _mm_mul_pd (_mm_mul_pd(d2a[i],   b[j]),   c[k]);
        d2_abc[1]   = _mm_mul_pd (_mm_mul_pd( da[i],  db[j]),   c[k]);
        d2_abc[2]   = _mm_mul_pd (_mm_mul_pd( da[i],   b[j]),  dc[k]);
        d2_abc[3]   = _mm_mul_pd (_mm_mul_pd(  a[i], d2b[j]),   c[k]);
        d2_abc[4]   = _mm_mul_pd (_mm_mul_pd(  a[i],  db[j]),  dc[k]);
        d2_abc[5]   = _mm_mul_pd (_mm_mul_pd(  a[i],   b[j]), d2c[k]);
        d3_abc[0] = _mm_mul_pd ( _mm_mul_pd(d3a[i],  b[j]),  c[k]);
        d3_abc[1] = _mm_mul_pd ( _mm_mul_pd(d2a[i], db[j]),  c[k]);
        d3_abc[2] = _mm_mul_pd ( _mm_mul_pd(d2a[i],  b[j]), dc[k]);
        d3_abc[3] = _mm_mul_pd ( _mm_mul_pd( da[i],d2b[j]),  c[k]);
        d3_abc[4] = _mm_mul_pd ( _mm_mul_pd( da[i], db[j]), dc[k]);
        d3_abc[5] = _mm_mul_pd ( _mm_mul_pd( da[i],  b[j]),d2c[k]);
        d3_abc[6] = _mm_mul_pd ( _mm_mul_pd(  a[i],d3b[j]),  c[k]);
        d3_abc[7] = _mm_mul_pd ( _mm_mul_pd(  a[i],d2b[j]), dc[k]);
        d3_abc[8] = _mm_mul_pd ( _mm_mul_pd(  a[i], db[j]),d2c[k]);
        d3_abc[9] = _mm_mul_pd ( _mm_mul_pd(  a[i],  b[j]),d3c[k]);
        __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<Nh; n++)
        {
          mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
          mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
          mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
          mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
          mhess[6*n+0]  = _mm_add_pd (mhess[6*n+0],  _mm_mul_pd (d2_abc[0], coefs[n]));
          mhess[6*n+1]  = _mm_add_pd (mhess[6*n+1],  _mm_mul_pd (d2_abc[1], coefs[n]));
          mhess[6*n+2]  = _mm_add_pd (mhess[6*n+2],  _mm_mul_pd (d2_abc[2], coefs[n]));
          mhess[6*n+3]  = _mm_add_pd (mhess[6*n+3],  _mm_mul_pd (d2_abc[3], coefs[n]));
          mhess[6*n+4]  = _mm_add_pd (mhess[6*n+4],  _mm_mul_pd (d2_abc[4], coefs[n]));
          mhess[6*n+5]  = _mm_add_pd (mhess[6*n+5],  _mm_mul_pd (d2_abc[5], coefs[n]));
          mgradhess[10*n+0] = _mm_add_pd (mgradhess[10*n+0], _mm_mul_pd(d3_abc[0], coefs[n]));
          mgradhess[10*n+1] = _mm_add_pd (mgradhess[10*n+1], _mm_mul_pd(d3_abc[1], coefs[n]));
          mgradhess[10*n+2] = _mm_add_pd (mgradhess[10*n+2], _mm_mul_pd(d3_abc[2], coefs[n]));
          mgradhess[10*n+3] = _mm_add_pd (mgradhess[10*n+3], _mm_mul_pd(d3_abc[3], coefs[n]));
          mgradhess[10*n+4] = _mm_add_pd (mgradhess[10*n+4], _mm_mul_pd(d3_abc[4], coefs[n]));
          mgradhess[10*n+5] = _mm_add_pd (mgradhess[10*n+5], _mm_mul_pd(d3_abc[5], coefs[n]));
          mgradhess[10*n+6] = _mm_add_pd (mgradhess[10*n+6], _mm_mul_pd(d3_abc[6], coefs[n]));
          mgradhess[10*n+7] = _mm_add_pd (mgradhess[10*n+7], _mm_mul_pd(d3_abc[7], coefs[n]));
          mgradhess[10*n+8] = _mm_add_pd (mgradhess[10*n+8], _mm_mul_pd(d3_abc[8], coefs[n]));
          mgradhess[10*n+9] = _mm_add_pd (mgradhess[10*n+9], _mm_mul_pd(d3_abc[9], coefs[n]));
        }
      }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<N/2; n++)
  {
    _mm_storeu_pd((double*)(vals+2*n),mvals[n]);
    _mm_storel_pd((double*)(grads+6*n+0),mgrads[3*n+0]);
    _mm_storeh_pd((double*)(grads+6*n+3),mgrads[3*n+0]);
    _mm_storel_pd((double*)(grads+6*n+1),mgrads[3*n+1]);
    _mm_storeh_pd((double*)(grads+6*n+4),mgrads[3*n+1]);
    _mm_storel_pd((double*)(grads+6*n+2),mgrads[3*n+2]);
    _mm_storeh_pd((double*)(grads+6*n+5),mgrads[3*n+2]);
    _mm_storel_pd((double*)(hess+18*n+0),  mhess [6*n+0]);
    _mm_storeh_pd((double*)(hess+18*n+9),  mhess [6*n+0]);
    _mm_storel_pd((double*)(hess+18*n+1),  mhess [6*n+1]);
    _mm_storeh_pd((double*)(hess+18*n+10), mhess [6*n+1]);
    _mm_storel_pd((double*)(hess+18*n+2),  mhess [6*n+2]);
    _mm_storeh_pd((double*)(hess+18*n+11), mhess [6*n+2]);
    _mm_storel_pd((double*)(hess+18*n+4),  mhess [6*n+3]);
    _mm_storeh_pd((double*)(hess+18*n+13), mhess [6*n+3]);
    _mm_storel_pd((double*)(hess+18*n+5),  mhess [6*n+4]);
    _mm_storeh_pd((double*)(hess+18*n+14), mhess [6*n+4]);
    _mm_storel_pd((double*)(hess+18*n+8),  mhess [6*n+5]);
    _mm_storeh_pd((double*)(hess+18*n+17), mhess [6*n+5]);
    _mm_storel_pd((double*)(gradhess+54*n+0),  mgradhess [10*n+0]);
    _mm_storeh_pd((double*)(gradhess+54*n+27), mgradhess [10*n+0]);
    _mm_storel_pd((double*)(gradhess+54*n+1),  mgradhess [10*n+1]);
    _mm_storeh_pd((double*)(gradhess+54*n+28), mgradhess [10*n+1]);
    _mm_storel_pd((double*)(gradhess+54*n+2),  mgradhess [10*n+2]);
    _mm_storeh_pd((double*)(gradhess+54*n+29), mgradhess [10*n+2]);
    _mm_storel_pd((double*)(gradhess+54*n+4),  mgradhess [10*n+3]);
    _mm_storeh_pd((double*)(gradhess+54*n+31), mgradhess [10*n+3]);
    _mm_storel_pd((double*)(gradhess+54*n+5),  mgradhess [10*n+4]);
    _mm_storeh_pd((double*)(gradhess+54*n+32), mgradhess [10*n+4]);
    _mm_storel_pd((double*)(gradhess+54*n+8),  mgradhess [10*n+5]);
    _mm_storeh_pd((double*)(gradhess+54*n+35), mgradhess [10*n+5]);
    _mm_storel_pd((double*)(gradhess+54*n+13), mgradhess [10*n+6]);
    _mm_storeh_pd((double*)(gradhess+54*n+40), mgradhess [10*n+6]);
    _mm_storel_pd((double*)(gradhess+54*n+14), mgradhess [10*n+7]);
    _mm_storeh_pd((double*)(gradhess+54*n+41), mgradhess [10*n+7]);
    _mm_storel_pd((double*)(gradhess+54*n+17), mgradhess [10*n+8]);
    _mm_storeh_pd((double*)(gradhess+54*n+44), mgradhess [10*n+8]);
    _mm_storel_pd((double*)(gradhess+54*n+26), mgradhess [10*n+9]);
    _mm_storeh_pd((double*)(gradhess+54*n+53), mgradhess [10*n+9]);
  }
  if (N&1)
  {
    _mm_storel_pd((double*)(vals+N-1),mvals[Nh-1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+0),mgrads[3*(Nh-1)+0]);
    _mm_storel_pd((double*)(grads+3*(N-1)+1),mgrads[3*(Nh-1)+1]);
    _mm_storel_pd((double*)(grads+3*(N-1)+2),mgrads[3*(Nh-1)+2]);
    _mm_storel_pd((double*)(hess+9*(N-1)+0),  mhess [6*(Nh-1)+0]);
    _mm_storel_pd((double*)(hess+9*(N-1)+1),  mhess [6*(Nh-1)+1]);
    _mm_storel_pd((double*)(hess+9*(N-1)+2),  mhess [6*(Nh-1)+2]);
    _mm_storel_pd((double*)(hess+9*(N-1)+4),  mhess [6*(Nh-1)+3]);
    _mm_storel_pd((double*)(hess+9*(N-1)+5),  mhess [6*(Nh-1)+4]);
    _mm_storel_pd((double*)(hess+9*(N-1)+8),  mhess [6*(Nh-1)+5]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+0),  mgradhess [10*(Nh-1)+0]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+1),  mgradhess [10*(Nh-1)+1]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+2),  mgradhess [10*(Nh-1)+2]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+4),  mgradhess [10*(Nh-1)+3]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+5),  mgradhess [10*(Nh-1)+4]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+8),  mgradhess [10*(Nh-1)+5]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+13), mgradhess [10*(Nh-1)+6]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+14), mgradhess [10*(Nh-1)+7]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+17), mgradhess [10*(Nh-1)+8]);
    _mm_storel_pd((double*)(gradhess+27*(Nh-1)+26), mgradhess [10*(Nh-1)+9]);
  }
  for (int n=0; n<N; n++)
  {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess[9*n+0]  *= dxInv*dxInv;
    hess[9*n+4]  *= dyInv*dyInv;
    hess[9*n+8]  *= dzInv*dzInv;
    hess[9*n+1]  *= dxInv*dyInv;
    hess[9*n+2]  *= dxInv*dzInv;
    hess[9*n+5]  *= dyInv*dzInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess[9*n+3] = hess[9*n+1];
    hess[9*n+6] = hess[9*n+2];
    hess[9*n+7] = hess[9*n+5];
    gradhess [27*n+0 ] *= dxInv*dxInv*dxInv;
    gradhess [27*n+1 ] *= dxInv*dxInv*dyInv;
    gradhess [27*n+2 ] *= dxInv*dxInv*dzInv;
    gradhess [27*n+4 ] *= dxInv*dyInv*dyInv;
    gradhess [27*n+5 ] *= dxInv*dyInv*dzInv;
    gradhess [27*n+8 ] *= dxInv*dzInv*dzInv;
    gradhess [27*n+13] *= dyInv*dyInv*dyInv;
    gradhess [27*n+14] *= dyInv*dyInv*dzInv;
    gradhess [27*n+17] *= dyInv*dzInv*dzInv;
    gradhess [27*n+26] *= dzInv*dzInv*dzInv;
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
