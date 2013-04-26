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

#ifndef BSPLINE_EVAL_SSE_Z_H
#define BSPLINE_EVAL_SSE_Z_H

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#include <math.h>

extern __m128d *restrict A_d;
// extern __m128d
//     A0_01,   A0_23,   A1_01,   A1_23,   A2_01,   A2_23,   A3_01,   A3_23,
//    dA0_01,  dA0_23,  dA1_01,  dA1_23,  dA2_01,  dA2_23,  dA3_01,  dA3_23,
//   d2A0_01, d2A0_23, d2A1_01, d2A1_23, d2A2_01, d2A2_23, d2A3_01, d2A3_23;



// This returns, pack in r, the two four-element dot products given
// by, r = [dot([a0,a1],[b0,b1], dot([a2,a3],[b2,b3]).  Specifically
// r_l = a0_l*b0_l + a0_h+b0_h + a1_l*b1_l + a1_h*b1_h
// r_h = a2_l*b2_l + a2_h+b2_h + a3_l*b1_l + a3_h*b1_h
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
inline void
eval_UBspline_1d_z (UBspline_1d_z * restrict spline,
                    double x, complex_double* restrict val)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  double tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;
  *val =
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_z_vg (UBspline_1d_z * restrict spline, double x,
                       complex_double* restrict val, complex_double* restrict grad)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  double tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;
  *val =
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv *
          (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
           coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
           coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
           coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_z_vgl (UBspline_1d_z * restrict spline, double x,
                        complex_double* restrict val, complex_double* restrict grad,
                        complex_double* restrict lapl)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  double tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;
  *val =
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv *
          (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
           coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
           coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
           coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
          (coefs[i+0]*(d2Ad[ 2]*tp[2] + d2Ad[ 3]*tp[3])+
           coefs[i+1]*(d2Ad[ 6]*tp[2] + d2Ad[ 7]*tp[3])+
           coefs[i+2]*(d2Ad[10]*tp[2] + d2Ad[11]*tp[3])+
           coefs[i+3]*(d2Ad[14]*tp[2] + d2Ad[15]*tp[3]));
}

inline void
eval_UBspline_1d_z_vgh (UBspline_1d_z * restrict spline, double x,
                        complex_double* restrict val, complex_double* restrict grad,
                        complex_double* restrict hess)
{
  eval_UBspline_1d_z_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_z (UBspline_2d_z * restrict spline,
                    double x, double y, complex_double* restrict val)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01, b01, a23, b23, bP01r, bP23r, bP01i, bP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
  tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
  tpx23 = _mm_set_pd (tx, 1.0);
  tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
  tpy23 = _mm_set_pd (ty, 1.0);
  // x-dependent vectors
  _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpx01, tpx23, tpx01, tpx23, a01);
  _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpx01, tpx23, tpx01, tpx23, a23);
  // y-dependent vectors
  _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpy01, tpy23, tpy01, tpy23, b01);
  _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpy01, tpy23, tpy01, tpy23, b23);
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_load_pd (P(0,0));
  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));
  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));
  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));
  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  tmp0 = _mm_load_pd (P(2,0));
  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));
  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));
  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));
  tmp1 = _mm_load_pd (P(3,3));
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
eval_UBspline_2d_z_vg (UBspline_2d_z * restrict spline,
                       double x, double y,
                       complex_double* restrict val, complex_double* restrict grad)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01, b01, da01, db01,
          a23, b23, da23, db23,
          bP01r, dbP01r,    bP23r, dbP23r,
          bP01i, dbP01i,    bP23i, dbP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  tmp0 = _mm_load_pd (P(0,0));
  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));
  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));
  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));
  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);
  tmp0 = _mm_load_pd (P(2,0));
  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));
  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));
  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));
  tmp1 = _mm_load_pd (P(3,3));
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
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_z_vgl (UBspline_2d_z * restrict spline,
                        double x, double y, complex_double* restrict val,
                        complex_double* restrict grad, complex_double* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01, b01, da01, db01, d2a01, d2b01,
          a23, b23, da23, db23, d2a23, d2b23,
          bP01r, dbP01r, d2bP01r,
          bP23r, dbP23r, d2bP23r,
          bP01i, dbP01i, d2bP01i,
          bP23i, dbP23i, d2bP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  tmp0 = _mm_load_pd (P(0,0));
  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));
  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));
  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));
  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01i);
  tmp0 = _mm_load_pd (P(2,0));
  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));
  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));
  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));
  tmp1 = _mm_load_pd (P(3,3));
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
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  d2x_r *= dxInv*dxInv;
  d2x_i *= dxInv*dxInv;
  d2y_r *= dyInv*dyInv;
  d2y_i *= dyInv*dyInv;
#ifdef __cplusplus
  *lapl = std::complex<double>(d2x_r + d2y_r, d2x_i + d2y_i);
#else
  *lapl = (d2x_r + d2y_r) + 1.0I*(d2x_i + d2y_i);
#endif
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_z_vgh (UBspline_2d_z * restrict spline,
                        double x, double y, complex_double* restrict val,
                        complex_double* restrict grad, complex_double* restrict hess)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23,
          a01, b01, da01, db01, d2a01, d2b01,
          a23, b23, da23, db23, d2a23, d2b23,
          bP01r, dbP01r, d2bP01r,
          bP23r, dbP23r, d2bP23r,
          bP01i, dbP01i, d2bP01i,
          bP23i, dbP23i, d2bP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  tmp0 = _mm_load_pd (P(0,0));
  tmp1 = _mm_load_pd (P(0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2));
  tmp1 = _mm_load_pd (P(0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0));
  tmp1 = _mm_load_pd (P(1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2));
  tmp1 = _mm_load_pd (P(1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,     b01,    b23,    b01,    b23,    bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,     b01,    b23,    b01,    b23,    bP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,    db01,   db23,   db01,   db23,   dbP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,    db01,   db23,   db01,   db23,   dbP01i);
  _MM_DDOT4_PD(r0, r1, r2, r3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01r);
  _MM_DDOT4_PD(i0, i1, i2, i3,   d2b01,  d2b23,  d2b01,  d2b23,  d2bP01i);
  tmp0 = _mm_load_pd (P(2,0));
  tmp1 = _mm_load_pd (P(2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2));
  tmp1 = _mm_load_pd (P(2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0));
  tmp1 = _mm_load_pd (P(3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2));
  tmp1 = _mm_load_pd (P(3,3));
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
  _MM_DOT4_PD (da01, da23, dbP01r, dbP23r, dhess[4]);
  _MM_DOT4_PD (da01, da23, dbP01i, dbP23i, dhess[5]);
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
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
eval_UBspline_3d_z (UBspline_3d_z * restrict spline,
                    double x, double y, double z,
                    complex_double* restrict val)
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
  int xs = spline->x_stride;
  int ys = spline->y_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01, b01, c01, cPr[8], dcPr[8], bcP01r, bcP23r,
          a23, b23, c23, cPi[8], dcPi[8], bcP01i, bcP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  tmp0 = _mm_load_pd (P(0,0,0));
  tmp1 = _mm_load_pd (P(0,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,0,2));
  tmp1 = _mm_load_pd (P(0,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,0));
  tmp1 = _mm_load_pd (P(0,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,2));
  tmp1 = _mm_load_pd (P(0,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  // 2nd eighth
  tmp0 = _mm_load_pd (P(0,2,0));
  tmp1 = _mm_load_pd (P(0,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2,2));
  tmp1 = _mm_load_pd (P(0,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,0));
  tmp1 = _mm_load_pd (P(0,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,2));
  tmp1 = _mm_load_pd (P(0,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  // 3rd eighth
  tmp0 = _mm_load_pd (P(1,0,0));
  tmp1 = _mm_load_pd (P(1,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0,2));
  tmp1 = _mm_load_pd (P(1,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,0));
  tmp1 = _mm_load_pd (P(1,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,2));
  tmp1 = _mm_load_pd (P(1,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  // 4th eighth
  tmp0 = _mm_load_pd (P(1,2,0));
  tmp1 = _mm_load_pd (P(1,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2,2));
  tmp1 = _mm_load_pd (P(1,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,0));
  tmp1 = _mm_load_pd (P(1,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,2));
  tmp1 = _mm_load_pd (P(1,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  // 5th eighth
  tmp0 = _mm_load_pd (P(2,0,0));
  tmp1 = _mm_load_pd (P(2,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,0,2));
  tmp1 = _mm_load_pd (P(2,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,0));
  tmp1 = _mm_load_pd (P(2,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,2));
  tmp1 = _mm_load_pd (P(2,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  // 6th eighth
  tmp0 = _mm_load_pd (P(2,2,0));
  tmp1 = _mm_load_pd (P(2,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2,2));
  tmp1 = _mm_load_pd (P(2,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,0));
  tmp1 = _mm_load_pd (P(2,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,2));
  tmp1 = _mm_load_pd (P(2,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  // 7th eighth
  tmp0 = _mm_load_pd (P(3,0,0));
  tmp1 = _mm_load_pd (P(3,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0,2));
  tmp1 = _mm_load_pd (P(3,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,0));
  tmp1 = _mm_load_pd (P(3,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,2));
  tmp1 = _mm_load_pd (P(3,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  // 8th eighth
  tmp0 = _mm_load_pd (P(3,2,0));
  tmp1 = _mm_load_pd (P(3,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2,2));
  tmp1 = _mm_load_pd (P(3,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,0));
  tmp1 = _mm_load_pd (P(3,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,2));
  tmp1 = _mm_load_pd (P(3,3,3));
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
eval_UBspline_3d_z_vg (UBspline_3d_z * restrict spline,
                       double x, double y, double z,
                       complex_double* restrict val, complex_double* restrict grad)
{
  _mm_prefetch((const char*) &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[14],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  int ys = spline->y_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01, b01, c01, da01, db01, dc01,
          a23, b23, c23, da23, db23, dc23,
          cPr[8], dcPr[8],
          cPi[8], dcPi[8],
          bcP01r, dbcP01r, bdcP01r,
          bcP23r, dbcP23r, bdcP23r,
          bcP01i, dbcP01i, bdcP01i,
          bcP23i, dbcP23i, bdcP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  tmp0 = _mm_load_pd (P(0,0,0));
  tmp1 = _mm_load_pd (P(0,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,0,2));
  tmp1 = _mm_load_pd (P(0,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,0));
  tmp1 = _mm_load_pd (P(0,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,2));
  tmp1 = _mm_load_pd (P(0,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  // 2nd eighth
  tmp0 = _mm_load_pd (P(0,2,0));
  tmp1 = _mm_load_pd (P(0,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2,2));
  tmp1 = _mm_load_pd (P(0,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,0));
  tmp1 = _mm_load_pd (P(0,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,2));
  tmp1 = _mm_load_pd (P(0,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  // 3rd eighth
  tmp0 = _mm_load_pd (P(1,0,0));
  tmp1 = _mm_load_pd (P(1,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0,2));
  tmp1 = _mm_load_pd (P(1,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,0));
  tmp1 = _mm_load_pd (P(1,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,2));
  tmp1 = _mm_load_pd (P(1,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);
  // 4th eighth
  tmp0 = _mm_load_pd (P(1,2,0));
  tmp1 = _mm_load_pd (P(1,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2,2));
  tmp1 = _mm_load_pd (P(1,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,0));
  tmp1 = _mm_load_pd (P(1,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,2));
  tmp1 = _mm_load_pd (P(1,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);
  // 5th eighth
  tmp0 = _mm_load_pd (P(2,0,0));
  tmp1 = _mm_load_pd (P(2,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,0,2));
  tmp1 = _mm_load_pd (P(2,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,0));
  tmp1 = _mm_load_pd (P(2,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,2));
  tmp1 = _mm_load_pd (P(2,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);
  // 6th eighth
  tmp0 = _mm_load_pd (P(2,2,0));
  tmp1 = _mm_load_pd (P(2,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2,2));
  tmp1 = _mm_load_pd (P(2,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,0));
  tmp1 = _mm_load_pd (P(2,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,2));
  tmp1 = _mm_load_pd (P(2,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);
  // 7th eighth
  tmp0 = _mm_load_pd (P(3,0,0));
  tmp1 = _mm_load_pd (P(3,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0,2));
  tmp1 = _mm_load_pd (P(3,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,0));
  tmp1 = _mm_load_pd (P(3,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,2));
  tmp1 = _mm_load_pd (P(3,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);
  // 8th eighth
  tmp0 = _mm_load_pd (P(3,2,0));
  tmp1 = _mm_load_pd (P(3,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2,2));
  tmp1 = _mm_load_pd (P(3,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,0));
  tmp1 = _mm_load_pd (P(3,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,2));
  tmp1 = _mm_load_pd (P(3,3,3));
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
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
#undef P
}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_z_vgl (UBspline_3d_z * restrict spline,
                        double x, double y, double z,
                        complex_double* restrict val, complex_double* restrict grad,
                        complex_double* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  int ys = spline->y_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,3), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
          a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
          cPr[8], dcPr[8], d2cPr[8],
          cPi[8], dcPi[8], d2cPi[8],
          bcP01r, dbcP01r, bdcP01r, d2bcP01r, bd2cP01r,
          bcP23r, dbcP23r, bdcP23r, d2bcP23r, bd2cP23r,
          bcP01i, dbcP01i, bdcP01i, d2bcP01i, bd2cP01i,
          bcP23i, dbcP23i, bdcP23i, d2bcP23i, bd2cP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  tmp0 = _mm_load_pd (P(0,0,0));
  tmp1 = _mm_load_pd (P(0,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,0,2));
  tmp1 = _mm_load_pd (P(0,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,0));
  tmp1 = _mm_load_pd (P(0,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,2));
  tmp1 = _mm_load_pd (P(0,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[0]);
  // 2nd eighth
  tmp0 = _mm_load_pd (P(0,2,0));
  tmp1 = _mm_load_pd (P(0,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2,2));
  tmp1 = _mm_load_pd (P(0,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,0));
  tmp1 = _mm_load_pd (P(0,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,2));
  tmp1 = _mm_load_pd (P(0,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[1]);
  // 3rd eighth
  tmp0 = _mm_load_pd (P(1,0,0));
  tmp1 = _mm_load_pd (P(1,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0,2));
  tmp1 = _mm_load_pd (P(1,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,0));
  tmp1 = _mm_load_pd (P(1,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,2));
  tmp1 = _mm_load_pd (P(1,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[2]);
  // 4th eighth
  tmp0 = _mm_load_pd (P(1,2,0));
  tmp1 = _mm_load_pd (P(1,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2,2));
  tmp1 = _mm_load_pd (P(1,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,0));
  tmp1 = _mm_load_pd (P(1,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,2));
  tmp1 = _mm_load_pd (P(1,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[3]);
  // 5th eighth
  tmp0 = _mm_load_pd (P(2,0,0));
  tmp1 = _mm_load_pd (P(2,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,0,2));
  tmp1 = _mm_load_pd (P(2,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,0));
  tmp1 = _mm_load_pd (P(2,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,2));
  tmp1 = _mm_load_pd (P(2,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[4]);
  // 6th eighth
  tmp0 = _mm_load_pd (P(2,2,0));
  tmp1 = _mm_load_pd (P(2,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2,2));
  tmp1 = _mm_load_pd (P(2,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,0));
  tmp1 = _mm_load_pd (P(2,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,2));
  tmp1 = _mm_load_pd (P(2,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[5]);
  // 7th eighth
  tmp0 = _mm_load_pd (P(3,0,0));
  tmp1 = _mm_load_pd (P(3,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0,2));
  tmp1 = _mm_load_pd (P(3,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,0));
  tmp1 = _mm_load_pd (P(3,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,2));
  tmp1 = _mm_load_pd (P(3,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[6]);
  // 8th eighth
  tmp0 = _mm_load_pd (P(3,2,0));
  tmp1 = _mm_load_pd (P(3,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2,2));
  tmp1 = _mm_load_pd (P(3,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,0));
  tmp1 = _mm_load_pd (P(3,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,2));
  tmp1 = _mm_load_pd (P(3,3,3));
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
  double sec_derivs[6];
  // Compute laplacian
  // d2x
  _MM_DOT4_PD (d2a01, d2a23, bcP01r, bcP23r, sec_derivs[0]);
  _MM_DOT4_PD (d2a01, d2a23, bcP01i, bcP23i, sec_derivs[1]);
  // d2y
  _MM_DOT4_PD (a01, a23, d2bcP01r, d2bcP23r, sec_derivs[2]);
  _MM_DOT4_PD (a01, a23, d2bcP01i, d2bcP23i, sec_derivs[3]);
  // d2z
  _MM_DOT4_PD (a01, a23, bd2cP01r, bd2cP23r, sec_derivs[4]);
  _MM_DOT4_PD (a01, a23, bd2cP01i, bd2cP23i, sec_derivs[5]);
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  sec_derivs[0] *= dxInv*dxInv;
  sec_derivs[1] *= dxInv*dxInv;
  sec_derivs[2] *= dyInv*dyInv;
  sec_derivs[3] *= dyInv*dyInv;
  sec_derivs[4] *= dzInv*dzInv;
  sec_derivs[5] *= dzInv*dzInv;
#ifdef __cplusplus
  *lapl = std::complex<double> (sec_derivs[0] + sec_derivs[2] + sec_derivs[4],
                                sec_derivs[1] + sec_derivs[3] + sec_derivs[5]);
#else
  *lapl = (sec_derivs[0] + sec_derivs[2] + sec_derivs[4]) +
          1.0I*(sec_derivs[1] + sec_derivs[3] + sec_derivs[5]);
#endif
#undef P
}



/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_z_vgh (UBspline_3d_z * restrict spline,
                        double x, double y, double z,
                        complex_double* restrict val, complex_double* restrict grad,
                        complex_double* restrict hess)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0);
  _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0);
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
  int xs = spline->x_stride;
  int ys = spline->y_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (const double*)(spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  _mm_prefetch ((const char*)P(0,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(0,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(1,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(2,3,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,0,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,1,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,2,3), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,0), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,1), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,2), _MM_HINT_T0);
  _mm_prefetch ((const char*)P(3,3,3), _MM_HINT_T0);
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
          a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
          a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
          cPr[8], dcPr[8], d2cPr[8],
          cPi[8], dcPi[8], d2cPi[8],
          bcP01r, dbcP01r, bdcP01r, d2bcP01r, dbdcP01r, bd2cP01r,
          bcP23r, dbcP23r, bdcP23r, d2bcP23r, dbdcP23r, bd2cP23r,
          bcP01i, dbcP01i, bdcP01i, d2bcP01i, dbdcP01i, bd2cP01i,
          bcP23i, dbcP23i, bdcP23i, d2bcP23i, dbdcP23i, bd2cP23i,
          tmp0, tmp1, tmp2, tmp3, r0, r1, r2, r3, i0, i1, i2, i3;
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
  // Compute cP, dcP, and d2cP products 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // Complex values are read in, then shuffled such that 4 registers
  // hold the read parts and 4 register hold the imaginary parts.
  // 1st eighth
  tmp0 = _mm_load_pd (P(0,0,0));
  tmp1 = _mm_load_pd (P(0,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,0,2));
  tmp1 = _mm_load_pd (P(0,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,0));
  tmp1 = _mm_load_pd (P(0,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,1,2));
  tmp1 = _mm_load_pd (P(0,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[0]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[0]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[0]);
  // 2nd eighth
  tmp0 = _mm_load_pd (P(0,2,0));
  tmp1 = _mm_load_pd (P(0,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,2,2));
  tmp1 = _mm_load_pd (P(0,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,0));
  tmp1 = _mm_load_pd (P(0,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(0,3,2));
  tmp1 = _mm_load_pd (P(0,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[1]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[1]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[1]);
  // 3rd eighth
  tmp0 = _mm_load_pd (P(1,0,0));
  tmp1 = _mm_load_pd (P(1,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,0,2));
  tmp1 = _mm_load_pd (P(1,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,0));
  tmp1 = _mm_load_pd (P(1,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,1,2));
  tmp1 = _mm_load_pd (P(1,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[2]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[2]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[2]);
  // 4th eighth
  tmp0 = _mm_load_pd (P(1,2,0));
  tmp1 = _mm_load_pd (P(1,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,2,2));
  tmp1 = _mm_load_pd (P(1,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,0));
  tmp1 = _mm_load_pd (P(1,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(1,3,2));
  tmp1 = _mm_load_pd (P(1,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[3]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[3]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[3]);
  // 5th eighth
  tmp0 = _mm_load_pd (P(2,0,0));
  tmp1 = _mm_load_pd (P(2,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,0,2));
  tmp1 = _mm_load_pd (P(2,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,0));
  tmp1 = _mm_load_pd (P(2,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,1,2));
  tmp1 = _mm_load_pd (P(2,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[4]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[4]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[4]);
  // 6th eighth
  tmp0 = _mm_load_pd (P(2,2,0));
  tmp1 = _mm_load_pd (P(2,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,2,2));
  tmp1 = _mm_load_pd (P(2,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,0));
  tmp1 = _mm_load_pd (P(2,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(2,3,2));
  tmp1 = _mm_load_pd (P(2,3,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[5]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[5]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[5]);
  // 7th eighth
  tmp0 = _mm_load_pd (P(3,0,0));
  tmp1 = _mm_load_pd (P(3,0,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,0,2));
  tmp1 = _mm_load_pd (P(3,0,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,0));
  tmp1 = _mm_load_pd (P(3,1,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,1,2));
  tmp1 = _mm_load_pd (P(3,1,3));
  r3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i3 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  _MM_DDOT4_PD(r0, r1, r2, r3,   c01,  c23,  c01,  c23,  cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,   c01,  c23,  c01,  c23,  cPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3,  dc01, dc23, dc01, dc23, dcPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3,  dc01, dc23, dc01, dc23, dcPi[6]);
  _MM_DDOT4_PD(r0, r1, r2, r3, d2c01,d2c23,d2c01,d2c23,d2cPr[6]);
  _MM_DDOT4_PD(i0, i1, i2, i3, d2c01,d2c23,d2c01,d2c23,d2cPi[6]);
  // 8th eighth
  tmp0 = _mm_load_pd (P(3,2,0));
  tmp1 = _mm_load_pd (P(3,2,1));
  r0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i0 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,2,2));
  tmp1 = _mm_load_pd (P(3,2,3));
  r1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i1 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,0));
  tmp1 = _mm_load_pd (P(3,3,1));
  r2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(0, 0));
  i2 = _mm_shuffle_pd (tmp0, tmp1, _MM_SHUFFLE2(1, 1));
  tmp0 = _mm_load_pd (P(3,3,2));
  tmp1 = _mm_load_pd (P(3,3,3));
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
  // Compute hessian
  // d2x
  _MM_DOT4_PD (d2a01, d2a23, bcP01r, bcP23r, ((double*)hess)[0]);
  _MM_DOT4_PD (d2a01, d2a23, bcP01i, bcP23i, ((double*)hess)[1]);
  // d2y
  _MM_DOT4_PD (a01, a23, d2bcP01r, d2bcP23r, ((double*)hess)[8]);
  _MM_DOT4_PD (a01, a23, d2bcP01i, d2bcP23i, ((double*)hess)[9]);
  // d2z
  _MM_DOT4_PD (a01, a23, bd2cP01r, bd2cP23r, ((double*)hess)[16]);
  _MM_DOT4_PD (a01, a23, bd2cP01i, bd2cP23i, ((double*)hess)[17]);
  // dx dy
  _MM_DOT4_PD (da01, da23, dbcP01r, dbcP23r, ((double*)hess)[2]);
  _MM_DOT4_PD (da01, da23, dbcP01i, dbcP23i, ((double*)hess)[3]);
  _MM_DOT4_PD (da01, da23, dbcP01r, dbcP23r, ((double*)hess)[6]);
  _MM_DOT4_PD (da01, da23, dbcP01i, dbcP23i, ((double*)hess)[7]);
  // dx dz
  _MM_DOT4_PD (da01, da23, bdcP01r, bdcP23r, ((double*)hess)[4]);
  _MM_DOT4_PD (da01, da23, bdcP01i, bdcP23i, ((double*)hess)[5]);
  _MM_DOT4_PD (da01, da23, bdcP01r, bdcP23r, ((double*)hess)[12]);
  _MM_DOT4_PD (da01, da23, bdcP01i, bdcP23i, ((double*)hess)[13]);
  // dy dz
  _MM_DOT4_PD (a01, a23, dbdcP01r, dbdcP23r, ((double*)hess)[10]);
  _MM_DOT4_PD (a01, a23, dbdcP01i, dbdcP23i, ((double*)hess)[11]);
  _MM_DOT4_PD (a01, a23, dbdcP01r, dbdcP23r, ((double*)hess)[14]);
  _MM_DOT4_PD (a01, a23, dbdcP01i, dbdcP23i, ((double*)hess)[15]);
  // Multiply gradients and hessians by appropriate grid inverses
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  hess[0] *= dxInv*dxInv;
  hess[4] *= dyInv*dyInv;
  hess[8] *= dzInv*dzInv;
  hess[1] *= dxInv*dyInv;
  hess[3] *= dxInv*dyInv;
  hess[2] *= dxInv*dzInv;
  hess[6] *= dxInv*dzInv;
  hess[5] *= dyInv*dzInv;
  hess[7] *= dyInv*dzInv;
  // Copy hessian elements into lower half of 3x3 matrix
  // hess[3] = hess[1];
  // hess[6] = hess[2];
  // hess[7] = hess[5];
#undef P
}

#endif
