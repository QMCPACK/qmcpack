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

#ifndef BSPLINE_EVAL_SSE_D_H
#define BSPLINE_EVAL_SSE_D_H

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
  #include <pmmintrin.h>
#endif
#include <math.h>

// extern __m128d 
//     A0_01,   A0_23,   A1_01,   A1_23,   A2_01,   A2_23,   A3_01,   A3_23, 
//    dA0_01,  dA0_23,  dA1_01,  dA1_23,  dA2_01,  dA2_23,  dA3_01,  dA3_23, 
//   d2A0_01, d2A0_23, d2A1_01, d2A1_23, d2A2_01, d2A2_23, d2A3_01, d2A3_23;

extern __m128d *restrict A_d;


extern double* restrict Ad;
extern double* restrict dAd;
extern double* restrict d2Ad;


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
/* NOTE:  SSE does not seem to speed things up in 1D.       */
/* Therefore, we simply copy the std routines.              */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_d (UBspline_1d_d * restrict spline, 
		    double x, double* restrict val)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  
  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_d_vg (UBspline_1d_d * restrict spline, double x, 
		     double* restrict val, double* restrict grad)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  
  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* restrict coefs = spline->coefs;

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
eval_UBspline_1d_d_vgl (UBspline_1d_d * restrict spline, double x, 
			double* restrict val, double* restrict grad,
			double* restrict lapl)
{
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;
  
  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* restrict coefs = spline->coefs;

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
eval_UBspline_1d_d_vgh (UBspline_1d_d * restrict spline, double x, 
			double* restrict val, double* restrict grad,
			double* restrict hess)
{
  eval_UBspline_1d_d_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_d (UBspline_2d_d * restrict spline, 
		    double x, double y, double* restrict val)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, 
    a01, b01, bP01, a23, b23, bP23,
    tmp0, tmp1, tmp2, tmp3;
  
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
  tmp0 = _mm_loadu_pd (P(0,0)); tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0)); tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  tmp0 = _mm_loadu_pd (P(2,0)); tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0)); tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
}


/* Value and gradient */
inline void
eval_UBspline_2d_d_vg (UBspline_2d_d * restrict spline, 
		       double x, double y, 
		       double* restrict val, double* restrict grad)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
  // Now compute the vectors:
  // tpx = [t_x^3 t_x^2 t_x 1]
  // tpy = [t_y^3 t_y^2 t_y 1]
  // tpz = [t_z^3 t_z^2 t_z 1]
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  __m128d tpx01, tpx23, tpy01, tpy23, 
    a01, b01, da01, db01, bP01, dbP01,  
    a23, b23, da23, db23, bP23, dbP23, 
    tmp0, tmp1, tmp2, tmp3;
  
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
  
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0)); tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0)); tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  tmp0 = _mm_loadu_pd (P(2,0)); tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0)); tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
#undef P
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_d_vgl (UBspline_2d_d * restrict spline, 
			double x, double y, double* restrict val, 
			double* restrict grad, double* restrict lapl)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
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
    bP01, dbP01, d2bP01, 
    bP23, dbP23, d2bP23,
    tmp0, tmp1, tmp2, tmp3;
  
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
  
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0)); tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0)); tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP01);
  tmp0 = _mm_loadu_pd (P(2,0)); tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0)); tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP23);
  
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
  // Compute laplacian
  double sec_derivs[2];
  _MM_DOT4_PD (d2a01, d2a23, bP01, bP23, sec_derivs[0]);
  _MM_DOT4_PD (a01, a23, d2bP01, d2bP23, sec_derivs[1]);
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  sec_derivs[0] *= dxInv * dxInv;
  sec_derivs[1] *= dyInv * dyInv;
  *lapl = sec_derivs[0] + sec_derivs[1];
#undef P
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_d_vgh (UBspline_2d_d * restrict spline, 
			double x, double y, double* restrict val, 
			double* restrict grad, double* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  int xs = spline->x_stride;
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j)))
//   _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
//   _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);

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
    bP01, dbP01, d2bP01, 
    bP23, dbP23, d2bP23,
    tmp0, tmp1, tmp2, tmp3;
  
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
  
  // Now compute bP, dbP, d2bP products
  tmp0 = _mm_loadu_pd (P(0,0)); tmp1 = _mm_loadu_pd(P(0,2));
  tmp2 = _mm_loadu_pd (P(1,0)); tmp3 = _mm_loadu_pd(P(1,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP01);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP01);
  tmp0 = _mm_loadu_pd (P(2,0)); tmp1 = _mm_loadu_pd(P(2,2));
  tmp2 = _mm_loadu_pd (P(3,0)); tmp3 = _mm_loadu_pd(P(3,2));
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,   b01,   b23,   b01,   b23,   bP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3,  db01,  db23,  db01,  db23,  dbP23);
  _MM_DDOT4_PD (tmp0, tmp1, tmp2, tmp3, d2b01, d2b23, d2b01, d2b23, d2bP23);
  
  // Compute value
  _MM_DOT4_PD (a01, a23, bP01, bP23, *val);
  // Compute gradient
  _MM_DOT4_PD (da01, da23, bP01, bP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbP01, dbP23, grad[1]);
  // Compute hessian
  _MM_DOT4_PD (d2a01, d2a23, bP01, bP23, hess[0]);
  _MM_DOT4_PD (a01, a23, d2bP01, d2bP23, hess[3]);
  _MM_DOT4_PD (da01, da23, dbP01, dbP23, hess[1]);
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  hess[0] *= dxInv * dxInv;
  hess[1] *= dxInv * dyInv;
  hess[3] *= dyInv * dyInv;
  hess[2] = hess[1];
#undef P
}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_d (UBspline_3d_d * restrict spline, 
		    double x, double y, double z,
		    double* restrict val)
{
  _mm_prefetch ((const char*)  &A_d[0],_MM_HINT_T0); _mm_prefetch ((const char*)  &A_d[1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[2],_MM_HINT_T0); _mm_prefetch ((const char*)  &A_d[3],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[4],_MM_HINT_T0); _mm_prefetch ((const char*)  &A_d[5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[6],_MM_HINT_T0); _mm_prefetch ((const char*)  &A_d[7],_MM_HINT_T0);  

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
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
    a01, b01, c01, a23, b23, c23,
    cP[8], dcP[8], d2cP[8], 
    bcP01, dbcP01, bdcP01, bcP23, dbcP23, bdcP23, 
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
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

  // Compute cP product 1/8 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st eighth
  tmp0    = _mm_loadu_pd (P(0,0,0));  tmp1    = _mm_loadu_pd (P(0,0,2));  
  tmp2    = _mm_loadu_pd (P(0,1,0));  tmp3    = _mm_loadu_pd (P(0,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  // 2nd eighth
  tmp0    = _mm_loadu_pd (P(0,2,0));  tmp1    = _mm_loadu_pd (P(0,2,2));
  tmp2    = _mm_loadu_pd (P(0,3,0));  tmp3    = _mm_loadu_pd (P(0,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  // 3rd eighth
  tmp0    = _mm_loadu_pd (P(1,0,0));  tmp1    = _mm_loadu_pd (P(1,0,2));
  tmp2    = _mm_loadu_pd (P(1,1,0));  tmp3    = _mm_loadu_pd (P(1,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  // 4th eighth
  tmp0    = _mm_loadu_pd (P(1,2,0));  tmp1    = _mm_loadu_pd (P(1,2,2));
  tmp2    = _mm_loadu_pd (P(1,3,0));  tmp3    = _mm_loadu_pd (P(1,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  // 5th eighth
  tmp0    = _mm_loadu_pd (P(2,0,0));  tmp1    = _mm_loadu_pd (P(2,0,2));
  tmp2    = _mm_loadu_pd (P(2,1,0));  tmp3    = _mm_loadu_pd (P(2,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  // 6th eighth
  tmp0    = _mm_loadu_pd (P(2,2,0));  tmp1    = _mm_loadu_pd (P(2,2,2));
  tmp2    = _mm_loadu_pd (P(2,3,0));  tmp3    = _mm_loadu_pd (P(2,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  // 7th eighth
  tmp0    = _mm_loadu_pd (P(3,0,0));  tmp1    = _mm_loadu_pd (P(3,0,2));
  tmp2    = _mm_loadu_pd (P(3,1,0));  tmp3    = _mm_loadu_pd (P(3,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  // 8th eighth
  tmp0    = _mm_loadu_pd (P(3,2,0));  tmp1    = _mm_loadu_pd (P(3,2,2));
  tmp2    = _mm_loadu_pd (P(3,3,0));  tmp3    = _mm_loadu_pd (P(3,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  
  // Now compute bcP products
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[0],   cP[1],   cP[2],   cP[3],   bcP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[4],   cP[5],   cP[6],   cP[7],   bcP23);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);

#undef P
}

/* Value and gradient */
inline void
eval_UBspline_3d_d_vg (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad)
{
  _mm_prefetch((const char*) &A_d[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[ 8],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[10],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[12],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);
  _mm_prefetch((const char*) &A_d[14],_MM_HINT_T0);  _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);

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
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
    a01, b01, c01, da01, db01, dc01, a23, b23, c23, da23, db23, dc23, 
    cP[8], dcP[8], d2cP[8], 
    bcP01, dbcP01, bdcP01, bcP23, dbcP23, bdcP23, 
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
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
  // 1st eighth
  tmp0    = _mm_loadu_pd (P(0,0,0));
  tmp1    = _mm_loadu_pd (P(0,0,2));
  tmp2    = _mm_loadu_pd (P(0,1,0));
  tmp3    = _mm_loadu_pd (P(0,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);

  // 2nd eighth
  tmp0    = _mm_loadu_pd (P(0,2,0));
  tmp1    = _mm_loadu_pd (P(0,2,2));
  tmp2    = _mm_loadu_pd (P(0,3,0));
  tmp3    = _mm_loadu_pd (P(0,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);


  // 3rd eighth
  tmp0    = _mm_loadu_pd (P(1,0,0));
  tmp1    = _mm_loadu_pd (P(1,0,2));
  tmp2    = _mm_loadu_pd (P(1,1,0));
  tmp3    = _mm_loadu_pd (P(1,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);

  // 4th eighth
  tmp0    = _mm_loadu_pd (P(1,2,0));
  tmp1    = _mm_loadu_pd (P(1,2,2));
  tmp2    = _mm_loadu_pd (P(1,3,0));
  tmp3    = _mm_loadu_pd (P(1,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);

  // 5th eighth
  tmp0    = _mm_loadu_pd (P(2,0,0));
  tmp1    = _mm_loadu_pd (P(2,0,2));
  tmp2    = _mm_loadu_pd (P(2,1,0));
  tmp3    = _mm_loadu_pd (P(2,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);

  // 6th eighth
  tmp0    = _mm_loadu_pd (P(2,2,0));
  tmp1    = _mm_loadu_pd (P(2,2,2));
  tmp2    = _mm_loadu_pd (P(2,3,0));
  tmp3    = _mm_loadu_pd (P(2,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);

  // 7th eighth
  tmp0    = _mm_loadu_pd (P(3,0,0));
  tmp1    = _mm_loadu_pd (P(3,0,2));
  tmp2    = _mm_loadu_pd (P(3,1,0));
  tmp3    = _mm_loadu_pd (P(3,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);

  // 8th eighth
  tmp0    = _mm_loadu_pd (P(3,2,0));
  tmp1    = _mm_loadu_pd (P(3,2,2));
  tmp2    = _mm_loadu_pd (P(3,3,0));
  tmp3    = _mm_loadu_pd (P(3,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[7]);

  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[0],   cP[1],   cP[2],   cP[3],   bcP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[4],   cP[5],   cP[6],   cP[7],   bcP23);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,   cP[0],   cP[1],   cP[2],   cP[3],  dbcP01);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,   cP[4],   cP[5],   cP[6],   cP[7],  dbcP23);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,  dcP[0],  dcP[1],  dcP[2],  dcP[3],  bdcP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,  dcP[4],  dcP[5],  dcP[6],  dcP[7],  bdcP23);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);

  // Compute gradient
  _MM_DOT4_PD (da01, da23, bcP01, bcP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbcP01, dbcP23, grad[1]);
  _MM_DOT4_PD (a01, a23, bdcP01, bdcP23, grad[2]);
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
eval_UBspline_3d_d_vgl (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad, double* restrict lapl)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

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
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
    a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
    a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
    cP[8], dcP[8], d2cP[8], 
    bcP01, dbcP01, bdcP01, d2bcP01, dbdcP01, bd2cP01,
    bcP23, dbcP23, bdcP23, d2bcP23, dbdcP23, bd2cP23,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
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
  // 1st eighth
  tmp0    = _mm_loadu_pd (P(0,0,0));
  tmp1    = _mm_loadu_pd (P(0,0,2));
  tmp2    = _mm_loadu_pd (P(0,1,0));
  tmp3    = _mm_loadu_pd (P(0,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[0]);

  // 2nd eighth
  tmp0    = _mm_loadu_pd (P(0,2,0));
  tmp1    = _mm_loadu_pd (P(0,2,2));
  tmp2    = _mm_loadu_pd (P(0,3,0));
  tmp3    = _mm_loadu_pd (P(0,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[1]);


  // 3rd eighth
  tmp0    = _mm_loadu_pd (P(1,0,0));
  tmp1    = _mm_loadu_pd (P(1,0,2));
  tmp2    = _mm_loadu_pd (P(1,1,0));
  tmp3    = _mm_loadu_pd (P(1,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[2]);

  // 4th eighth
  tmp0    = _mm_loadu_pd (P(1,2,0));
  tmp1    = _mm_loadu_pd (P(1,2,2));
  tmp2    = _mm_loadu_pd (P(1,3,0));
  tmp3    = _mm_loadu_pd (P(1,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[3]);

  // 5th eighth
  tmp0    = _mm_loadu_pd (P(2,0,0));
  tmp1    = _mm_loadu_pd (P(2,0,2));
  tmp2    = _mm_loadu_pd (P(2,1,0));
  tmp3    = _mm_loadu_pd (P(2,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[4]);

  // 6th eighth
  tmp0    = _mm_loadu_pd (P(2,2,0));
  tmp1    = _mm_loadu_pd (P(2,2,2));
  tmp2    = _mm_loadu_pd (P(2,3,0));
  tmp3    = _mm_loadu_pd (P(2,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[5]);

  // 7th eighth
  tmp0    = _mm_loadu_pd (P(3,0,0));
  tmp1    = _mm_loadu_pd (P(3,0,2));
  tmp2    = _mm_loadu_pd (P(3,1,0));
  tmp3    = _mm_loadu_pd (P(3,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[6]);

  // 8th eighth
  tmp0    = _mm_loadu_pd (P(3,2,0));
  tmp1    = _mm_loadu_pd (P(3,2,2));
  tmp2    = _mm_loadu_pd (P(3,3,0));
  tmp3    = _mm_loadu_pd (P(3,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[7]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[7]);

  
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[0],   cP[1],   cP[2],   cP[3],   bcP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,   cP[4],   cP[5],   cP[6],   cP[7],   bcP23);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,   cP[0],   cP[1],   cP[2],   cP[3],  dbcP01);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,   cP[4],   cP[5],   cP[6],   cP[7],  dbcP23);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,  dcP[0],  dcP[1],  dcP[2],  dcP[3],  bdcP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23,  dcP[4],  dcP[5],  dcP[6],  dcP[7],  bdcP23);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23,   cP[0],   cP[1],   cP[2],   cP[3], d2bcP01);
  _MM_DDOT4_PD (d2b01, d2b23, d2b01, d2b23,   cP[4],   cP[5],   cP[6],   cP[7], d2bcP23);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23, d2cP[0], d2cP[1], d2cP[2], d2cP[3], bd2cP01);
  _MM_DDOT4_PD (  b01,   b23,   b01,   b23, d2cP[4], d2cP[5], d2cP[6], d2cP[7], bd2cP23);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,  dcP[0],  dcP[1],  dcP[2],  dcP[3], dbdcP01);
  _MM_DDOT4_PD ( db01,  db23,  db01,  db23,  dcP[4],  dcP[5],  dcP[6],  dcP[7], dbdcP23);

  // Compute value
  _MM_DOT4_PD (a01, a23, bcP01, bcP23, *val);

  // Compute gradient
  _MM_DOT4_PD (da01, da23, bcP01, bcP23, grad[0]);
  _MM_DOT4_PD (a01, a23, dbcP01, dbcP23, grad[1]);
  _MM_DOT4_PD (a01, a23, bdcP01, bdcP23, grad[2]);
  // Compute laplacian
  double lx, ly, lz;
  _MM_DOT4_PD (d2a01, d2a23, bcP01, bcP23, lx);  // d2x
  _MM_DOT4_PD (a01, a23, d2bcP01, d2bcP23, ly);  // d2y
  _MM_DOT4_PD (a01, a23, bd2cP01, bd2cP23, lz);  // d2z
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  lx      *= dxInv*dxInv;
  ly      *= dyInv*dyInv;
  lz      *= dzInv*dzInv;
  *lapl = lx + ly + lz;
#undef P
}



/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_d_vgh (UBspline_3d_d * restrict spline, 
			double x, double y, double z,
			double* restrict val, double* restrict grad, 
			double* restrict hess)
{
  _mm_prefetch ((const char*)  &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
  _mm_prefetch ((const char*)  &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

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
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  int xs = spline->x_stride;
  int ys = spline->y_stride;

  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j,k) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz+k))
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
    a01, b01, c01, da01, db01, dc01, d2a01, d2b01, d2c01,
    a23, b23, c23, da23, db23, dc23, d2a23, d2b23, d2c23,
    cP[8], dcP[8], d2cP[8], 
    bcP01, dbcP01, bdcP01, d2bcP01, dbdcP01, bd2cP01,
    bcP23, dbcP23, bdcP23, d2bcP23, dbdcP23, bd2cP23,
    tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
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
  // 1st eighth
  tmp0    = _mm_loadu_pd (P(0,0,0));
  tmp1    = _mm_loadu_pd (P(0,0,2));
  tmp2    = _mm_loadu_pd (P(0,1,0));
  tmp3    = _mm_loadu_pd (P(0,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[0]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[0]);

  // 2nd eighth
  tmp0    = _mm_loadu_pd (P(0,2,0));
  tmp1    = _mm_loadu_pd (P(0,2,2));
  tmp2    = _mm_loadu_pd (P(0,3,0));
  tmp3    = _mm_loadu_pd (P(0,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[1]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[1]);


  // 3rd eighth
  tmp0    = _mm_loadu_pd (P(1,0,0));
  tmp1    = _mm_loadu_pd (P(1,0,2));
  tmp2    = _mm_loadu_pd (P(1,1,0));
  tmp3    = _mm_loadu_pd (P(1,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[2]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[2]);

  // 4th eighth
  tmp0    = _mm_loadu_pd (P(1,2,0));
  tmp1    = _mm_loadu_pd (P(1,2,2));
  tmp2    = _mm_loadu_pd (P(1,3,0));
  tmp3    = _mm_loadu_pd (P(1,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[3]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[3]);

  // 5th eighth
  tmp0    = _mm_loadu_pd (P(2,0,0));
  tmp1    = _mm_loadu_pd (P(2,0,2));
  tmp2    = _mm_loadu_pd (P(2,1,0));
  tmp3    = _mm_loadu_pd (P(2,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[4]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[4]);

  // 6th eighth
  tmp0    = _mm_loadu_pd (P(2,2,0));
  tmp1    = _mm_loadu_pd (P(2,2,2));
  tmp2    = _mm_loadu_pd (P(2,3,0));
  tmp3    = _mm_loadu_pd (P(2,3,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[5]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[5]);

  // 7th eighth
  tmp0    = _mm_loadu_pd (P(3,0,0));
  tmp1    = _mm_loadu_pd (P(3,0,2));
  tmp2    = _mm_loadu_pd (P(3,1,0));
  tmp3    = _mm_loadu_pd (P(3,1,2));
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,  c01,  c23,  c01,  c23,  cP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3, dc01, dc23, dc01, dc23, dcP[6]);
  _MM_DDOT4_PD(tmp0,tmp1,tmp2,tmp3,d2c01,d2c23,d2c01,d2c23,d2cP[6]);

  // 8th eighth
  tmp0    = _mm_loadu_pd (P(3,2,0));
  tmp1    = _mm_loadu_pd (P(3,2,2));
  tmp2    = _mm_loadu_pd (P(3,3,0));
  tmp3    = _mm_loadu_pd (P(3,3,2));
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
  hess[2] *= dxInv*dzInv;
  hess[5] *= dyInv*dzInv;
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];

#undef P


}


//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[0], tpx01), _mm_mul_pd (A_d[1], tpx23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[2], tpx01), _mm_mul_pd (A_d[3], tpx23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[4], tpx01), _mm_mul_pd (A_d[5], tpx23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[6], tpx01), _mm_mul_pd (A_d[7], tpx23));
//   a01  = _mm_hadd_pd(tmp0, tmp1);
//   a23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[8], tpx01), _mm_mul_pd (A_d[9], tpx23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[10], tpx01), _mm_mul_pd (A_d[11], tpx23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[12], tpx01), _mm_mul_pd (A_d[13], tpx23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[14], tpx01), _mm_mul_pd (A_d[15], tpx23));
//   da01  = _mm_hadd_pd(tmp0, tmp1);
//   da23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[16], tpx01), _mm_mul_pd (A_d[17], tpx23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[18], tpx01), _mm_mul_pd (A_d[19], tpx23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[20], tpx01), _mm_mul_pd (A_d[21], tpx23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[22], tpx01), _mm_mul_pd (A_d[23], tpx23));
//   d2a01  = _mm_hadd_pd(tmp0, tmp1);
//   d2a23  = _mm_hadd_pd(tmp2, tmp3);


//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[0], tpy01), _mm_mul_pd (A_d[1], tpy23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[2], tpy01), _mm_mul_pd (A_d[3], tpy23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[4], tpy01), _mm_mul_pd (A_d[5], tpy23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[6], tpy01), _mm_mul_pd (A_d[7], tpy23));
//   b01  = _mm_hadd_pd(tmp0, tmp1);
//   b23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[8], tpy01), _mm_mul_pd (A_d[9], tpy23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[10], tpy01), _mm_mul_pd (A_d[11], tpy23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[12], tpy01), _mm_mul_pd (A_d[13], tpy23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[14], tpy01), _mm_mul_pd (A_d[15], tpy23));
//   db01  = _mm_hadd_pd(tmp0, tmp1);
//   db23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[16], tpy01), _mm_mul_pd (A_d[17], tpz23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[18], tpy01), _mm_mul_pd (A_d[19], tpz23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[20], tpy01), _mm_mul_pd (A_d[21], tpz23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[22], tpy01), _mm_mul_pd (A_d[23], tpz23));
//   d2b01  = _mm_hadd_pd(tmp0, tmp1);
//   d2b23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[0], tpz01), _mm_mul_pd (A_d[1], tpz23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[2], tpz01), _mm_mul_pd (A_d[3], tpz23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[4], tpz01), _mm_mul_pd (A_d[5], tpz23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[6], tpz01), _mm_mul_pd (A_d[7], tpz23));
//   c01  = _mm_hadd_pd(tmp0, tmp1);
//   c23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[8], tpz01), _mm_mul_pd (A_d[9], tpz23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[10], tpz01), _mm_mul_pd (A_d[11], tpz23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[12], tpz01), _mm_mul_pd (A_d[13], tpz23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[14], tpz01), _mm_mul_pd (A_d[15], tpz23));
//   dc01  = _mm_hadd_pd(tmp0, tmp1);
//   dc23  = _mm_hadd_pd(tmp2, tmp3);

//   tmp0 = _mm_hadd_pd(_mm_mul_pd (A_d[16], tpz01), _mm_mul_pd (A_d[17], tpz23));
//   tmp1 = _mm_hadd_pd(_mm_mul_pd (A_d[18], tpz01), _mm_mul_pd (A_d[19], tpz23));
//   tmp2 = _mm_hadd_pd(_mm_mul_pd (A_d[20], tpz01), _mm_mul_pd (A_d[21], tpz23));
//   tmp3 = _mm_hadd_pd(_mm_mul_pd (A_d[22], tpz01), _mm_mul_pd (A_d[23], tpz23));
//   d2c01  = _mm_hadd_pd(tmp0, tmp1);
//   d2c23  = _mm_hadd_pd(tmp2, tmp3);
 
 


#endif
