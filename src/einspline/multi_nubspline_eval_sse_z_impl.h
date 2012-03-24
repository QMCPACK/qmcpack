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

#ifndef MULTI_BSPLINE_EVAL_SSE_Z_IMPL_H
#define MULTI_BSPLINE_EVAL_SSE_Z_IMPL_H

#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
#include <math.h>
#include "bspline_base.h"
#include "multi_nubspline_structs.h"

extern __m128d *restrict A_d;
extern double *restrict Ad, *restrict dAd, *restrict d2Ad;

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
/* 1D double-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_NUBspline_1d_z (multi_NUBspline_1d_z *spline,
			  double x,
			  complex_double* restrict vals)
{
  double a[4];

  int ix = get_NUBasis_funcs_d (spline->x_basis, x, a);
  int xs = spline->x_stride;

  complex_double* restrict coefs0 = spline->coefs +(ix+0)*xs;
  complex_double* restrict coefs1 = spline->coefs +(ix+1)*xs;
  complex_double* restrict coefs2 = spline->coefs +(ix+2)*xs;
  complex_double* restrict coefs3 = spline->coefs +(ix+3)*xs;
  for (int n=0; n<spline->num_splines; n++) 
    vals[n] = (a[0]*coefs0[n] + a[1]*coefs1[n] + 
	       a[2]*coefs2[n] + a[3]*coefs3[n]);
}



void
eval_multi_NUBspline_1d_z_vg (multi_NUBspline_1d_z *spline,
			     double x,
			     complex_double* restrict vals,
			     complex_double* restrict grads)
{
  double a[4], da[4];
  int ix = get_NUBasis_dfuncs_d (spline->x_basis, x, a, da);
  int xs = spline->x_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }

  for (int i=0; i<4; i++) { 
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++) {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
    }
  }
}


void
eval_multi_NUBspline_1d_z_vgl (multi_NUBspline_1d_z *spline,
			       double x,
			       complex_double* restrict vals,
			       complex_double* restrict grads,
			       complex_double* restrict lapl)	  
{
  double a[4], da[4], d2a[4];
  int ix = get_NUBasis_d2funcs_d (spline->x_basis, x, a, da, d2a);
  int xs = spline->x_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }

  for (int i=0; i<4; i++) {      
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++) {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
      lapl[n]  += d2a[i] * coefs[n];
    }
  }
}


void
eval_multi_NUBspline_1d_z_vgh (multi_NUBspline_1d_z *spline,
			      double x,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict hess)
{
  eval_multi_NUBspline_1d_z_vgl (spline, x, vals, grads, hess);
}



// /************************************************************/
// /* 2D double-precision, complex evaulation functions        */
// /************************************************************/
// void
// eval_multi_NUBspline_2d_z (multi_NUBspline_2d_z *spline,
// 			   double x, double y,
// 			   complex_double* restrict vals)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   double ipartx, iparty, tx, ty;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23,
//     a01  ,   b01,   a23,    b23;  

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);

//   // Zero-out values
//   __m128d mvals[N];
//   for (int n=0; n<N; n++) 
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
  
//   __m128d a[4], b[4];
//   a[0]=_mm_unpacklo_pd(a01,a01);
//   a[1]=_mm_unpackhi_pd(a01,a01);
//   a[2]=_mm_unpacklo_pd(a23,a23);
//   a[3]=_mm_unpackhi_pd(a23,a23);
				
//   b[0]=_mm_unpacklo_pd(b01,b01);
//   b[1]=_mm_unpackhi_pd(b01,b01);
//   b[2]=_mm_unpacklo_pd(b23,b23);
//   b[3]=_mm_unpackhi_pd(b23,b23);
				 				   				  
//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       __m128d ab              = _mm_mul_pd(  a[i],  b[j]);
//       __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);

//       for (int n=0; n<N; n++) 
// 	mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
//     }
  
//   for (int n=0; n<N; n++) 
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
// }



// void
// eval_multi_NUBspline_2d_z_vg (multi_NUBspline_2d_z *spline,
// 			     double x, double y,
// 			     complex_double* restrict vals,
// 			     complex_double* restrict grads)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   double ipartx, iparty, tx, ty;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23,
//     a01  ,   b01,   a23,    b23,  
//     da01 ,  db01,  da23,   db23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);

//   // Zero-out values
//   __m128d mvals[N], mgrads[2*N];
//   for (int n=0; n<N; n++) {
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
//     for (int i=0; i<2; i++) 
//       mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i],mgrads[2*n+i]);
//   }   
  
//   __m128d a[4], b[4], da[4], db[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01);
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01);
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23);
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23);
				 				  
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01);
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01);
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23);
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23);
				 				   				  
//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       __m128d ab, d_ab[2];
	
// 	ab         = _mm_mul_pd(  a[i],  b[j]);
// 	d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
// 	d_ab[1]    = _mm_mul_pd(  a[i], db[j]);

// 	__m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);

// 	for (int n=0; n<N; n++) {
// 	  mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
// 	  mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
// 	  mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
// 	}
//       }
  
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   complex_double lapl2[2*N];
//   for (int n=0; n<N; n++) {
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
//     _mm_storeu_pd((double*)(grads+2*n+0), mgrads[2*n+0]);
//     _mm_storeu_pd((double*)(grads+2*n+1), mgrads[2*n+1]);
//   }
//   for (int n=0; n<N; n++) {
//     grads[2*n+0] *= dxInv;
//     grads[2*n+1] *= dyInv;
//   }
// }



// void
// eval_multi_NUBspline_2d_z_vgl (multi_NUBspline_2d_z *spline,
// 			      double x, double y,
// 			      complex_double* restrict vals,
// 			      complex_double* restrict grads,
// 			      complex_double* restrict lapl)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   double ipartx, iparty, tx, ty;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23,
//     a01  ,   b01,   a23,    b23,  
//     da01 ,  db01,  da23,   db23,  
//     d2a01, d2b01, d2a23,  d2b23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);

//   // Zero-out values
//   __m128d mvals[N], mgrads[2*N], mlapl[2*N];
//   for (int n=0; n<N; n++) {
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
//     for (int i=0; i<2; i++) {
//       mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i],mgrads[2*n+i]);
//       mlapl [2*n+i]  = _mm_sub_pd (mlapl[2*n+i], mlapl[2*n+i]);
//     }
//   }   
  
//   __m128d a[4], b[4], da[4], db[4], d2a[4], d2b[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01); d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01); d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23); d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23); d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
				 				   				  
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01); d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01); d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23); d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23); d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
				 				   				  
//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       __m128d ab, d_ab[2], d2_ab[2];
	
// 	ab         = _mm_mul_pd(  a[i],  b[j]);
// 	d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
// 	d_ab[1]    = _mm_mul_pd(  a[i], db[j]);
// 	d2_ab[0]   = _mm_mul_pd(d2a[i],   b[j]);
// 	d2_ab[1]   = _mm_mul_pd(  a[i], d2b[j]);

// 	__m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);

// 	for (int n=0; n<N; n++) {
// 	  mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
// 	  mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
// 	  mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
// 	  mlapl [2*n+0] = _mm_add_pd (mlapl [2*n+0], _mm_mul_pd (d2_ab[0], coefs[n]));
// 	  mlapl [2*n+1] = _mm_add_pd (mlapl [2*n+1], _mm_mul_pd (d2_ab[1], coefs[n]));
// 	}
//       }
  
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   complex_double lapl2[2*N];
//   for (int n=0; n<N; n++) {
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
//     _mm_storeu_pd((double*)(grads+2*n+0), mgrads[2*n+0]);
//     _mm_storeu_pd((double*)(grads+2*n+1), mgrads[2*n+1]);
//     _mm_storeu_pd((double*)(lapl2+2*n+0), mlapl [2*n+0]);
//     _mm_storeu_pd((double*)(lapl2+2*n+1), mlapl [2*n+1]);
//   }
//   for (int n=0; n<N; n++) {
//     grads[2*n+0] *= dxInv;
//     grads[2*n+1] *= dyInv;
//     lapl2[2*n+0]  *= dxInv*dxInv;
//     lapl2[2*n+1]  *= dyInv*dyInv;
//     lapl[n] = lapl2[2*n+0] + lapl2[2*n+1];
//   }
// }




// void
// eval_multi_NUBspline_2d_z_vgh (multi_NUBspline_2d_z *spline,
// 			      double x, double y,
// 			      complex_double* restrict vals,
// 			      complex_double* restrict grads,
// 			      complex_double* restrict hess)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   double ipartx, iparty, tx, ty;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23,
//     a01  ,   b01,   a23,    b23,  
//     da01 ,  db01,  da23,   db23,  
//     d2a01, d2b01, d2a23,  d2b23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);

//   // Zero-out values
//   __m128d mvals[N], mgrads[2*N], mhess[3*N];
//   for (int n=0; n<N; n++) {
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
//     for (int i=0; i<2; i++)
//       mgrads[2*n+i] = _mm_sub_pd (mgrads[2*n+i],mgrads[2*n+i]);
//     for (int i=0; i<3; i++)
//       mhess[3*n+i]  = _mm_sub_pd (mhess[3*n+i], mhess[3*n+i]);
//   }   

//   __m128d a[4], b[4], da[4], db[4], d2a[4], d2b[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01); d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01); d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23); d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23); d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
				 				   				  
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01); d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01); d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23); d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23); d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
				 				   				  
//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       __m128d ab, d_ab[2], d2_ab[3];
	
// 	ab         = _mm_mul_pd(  a[i],  b[j]);
// 	d_ab[0]    = _mm_mul_pd(da[i],   b[j]);
// 	d_ab[1]    = _mm_mul_pd(  a[i], db[j]);
// 	d2_ab[0]   = _mm_mul_pd(d2a[i],   b[j]);
// 	d2_ab[1]   = _mm_mul_pd( da[i],  db[j]);
// 	d2_ab[2]   = _mm_mul_pd(  a[i], d2b[j]);

// 	__m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);

// 	for (int n=0; n<N; n++) {
// 	  mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   ab   , coefs[n]));
// 	  mgrads[2*n+0] = _mm_add_pd (mgrads[2*n+0], _mm_mul_pd ( d_ab[0], coefs[n]));
// 	  mgrads[2*n+1] = _mm_add_pd (mgrads[2*n+1], _mm_mul_pd ( d_ab[1], coefs[n]));
// 	  mhess[3*n+0]  = _mm_add_pd (mhess[3*n+0],  _mm_mul_pd (d2_ab[0], coefs[n]));
// 	  mhess[3*n+1]  = _mm_add_pd (mhess[3*n+1],  _mm_mul_pd (d2_ab[1], coefs[n]));
// 	  mhess[3*n+2]  = _mm_add_pd (mhess[3*n+2],  _mm_mul_pd (d2_ab[2], coefs[n]));
// 	}
//       }
  
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
  
//   for (int n=0; n<N; n++) {
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
//     _mm_storeu_pd((double*)(grads+2*n+0),mgrads[2*n+0]);
//     _mm_storeu_pd((double*)(grads+2*n+1),mgrads[2*n+1]);
//     _mm_storeu_pd((double*)(hess+4*n+0), mhess [3*n+0]);
//     _mm_storeu_pd((double*)(hess+4*n+1), mhess [3*n+1]);
//     _mm_storeu_pd((double*)(hess+4*n+3), mhess [3*n+2]);
//   }
//   for (int n=0; n<N; n++) {
//     grads[2*n+0] *= dxInv;
//     grads[2*n+1] *= dyInv;
//     hess[4*n+0]  *= dxInv*dxInv;
//     hess[4*n+1]  *= dxInv*dyInv;
//     hess[4*n+3]  *= dyInv*dyInv;
//     // Copy hessian elements into lower half of 3x3 matrix
//     hess[4*n+2] = hess[4*n+1];
//   }
// }


// /************************************************************/
// /* 3D double-precision, complex evaulation functions        */
// /************************************************************/
// void
// eval_multi_NUBspline_3d_z (multi_NUBspline_3d_z *spline,
// 			  double x, double y, double z,
// 			  complex_double* restrict vals)
// {
//   _mm_prefetch ((const char*) &A_d[0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[7],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   z -= spline->z_grid.start;
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   double uz = z*spline->z_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
//   double ipartx, iparty, ipartz, tx, ty, tz;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
//   tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int zs = spline->z_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
//     a01, b01, c01, a23, b23, c23,  
//     tmp0, tmp1, r0, r1, i0, i1, val_r, val_i;
  
//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);
//   tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
//   tpz23 = _mm_set_pd (tz, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpx01, tpx23, tpx01, tpx23,   a23);
//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpy01, tpy23, tpy01, tpy23,   b23);
//   // z-dependent vectors
//   _MM_DDOT4_PD (A_d[0], A_d[1], A_d[2], A_d[3], tpz01, tpz23, tpz01, tpz23,   c01);
//   _MM_DDOT4_PD (A_d[4], A_d[5], A_d[6], A_d[7], tpz01, tpz23, tpz01, tpz23,   c23);

//   // Zero-out values
//   __m128d mvals[N];
//   for (int n=0; n<N; n++) 
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);

//   __m128d a[4], b[4], c[4];
//   a[0]   = _mm_unpacklo_pd(a01,a01);
//   a[1]   = _mm_unpackhi_pd(a01,a01);
//   a[2]   = _mm_unpacklo_pd(a23,a23);
//   a[3]   = _mm_unpackhi_pd(a23,a23);

//   b[0]   = _mm_unpacklo_pd(b01,b01);
//   b[1]   = _mm_unpackhi_pd(b01,b01);
//   b[2]   = _mm_unpacklo_pd(b23,b23);
//   b[3]   = _mm_unpackhi_pd(b23,b23);

//   c[0]   = _mm_unpacklo_pd(c01,c01);
//   c[1]   = _mm_unpackhi_pd(c01,c01);
//   c[2]   = _mm_unpacklo_pd(c23,c23);
//   c[3]   = _mm_unpackhi_pd(c23,c23);

// #ifdef USE_PREFETCH
//   const int offset = PREFETCH_AHEAD;
// #else
//   const int offset = 0;
// #endif
//   int Nstop = N - offset;
//   if (Nstop & 1) 
//     Nstop--;
//   if (Nstop < 0)
//     Nstop = 0;

//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
      
// #ifdef USE_PREFETCH
//       __m128d abc[4];
//       for (int k=0; k<4; k++) 
// 	abc[k] = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
//       __m128d* restrict coefs0 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+0)*zs);
//       __m128d* restrict coefs1 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+1)*zs);
//       __m128d* restrict coefs2 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+2)*zs);
//       __m128d* restrict coefs3 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+3)*zs);
      
//       for (int n=0; n<Nstop; n+=2) { 
// 	//	_mm_prefetch ((const char*)&(mvals[n+offset]), _MM_HINT_NTA);
// 	_mm_prefetch ((const char*)&(coefs0[n+offset]), _MM_HINT_NTA);
// 	mvals[n+0] = _mm_add_pd(mvals[n+0], _mm_mul_pd (abc[0], coefs0[n+0]));
// 	mvals[n+0] = _mm_add_pd(mvals[n+0], _mm_mul_pd (abc[1], coefs1[n+0]));
// 	_mm_prefetch ((const char*)&(coefs1[n+offset]), _MM_HINT_NTA);
// 	mvals[n+0] = _mm_add_pd(mvals[n+0], _mm_mul_pd (abc[2], coefs2[n+0]));
// 	mvals[n+0] = _mm_add_pd(mvals[n+0], _mm_mul_pd (abc[3], coefs3[n+0]));
// 	_mm_prefetch ((const char*)&(coefs2[n+offset]), _MM_HINT_NTA);
// 	mvals[n+1] = _mm_add_pd(mvals[n+1], _mm_mul_pd (abc[0], coefs0[n+1]));
// 	mvals[n+1] = _mm_add_pd(mvals[n+1], _mm_mul_pd (abc[1], coefs1[n+1]));
// 	_mm_prefetch ((const char*)&(coefs3[n+offset]), _MM_HINT_NTA);
// 	mvals[n+1] = _mm_add_pd(mvals[n+1], _mm_mul_pd (abc[2], coefs2[n+1]));
// 	mvals[n+1] = _mm_add_pd(mvals[n+1], _mm_mul_pd (abc[3], coefs3[n+1]));
//       }
//       for (int n=Nstop; n<N; n++) {
// 	mvals[n] = _mm_add_pd(mvals[n], _mm_mul_pd (abc[0], coefs0[n]));
// 	mvals[n] = _mm_add_pd(mvals[n], _mm_mul_pd (abc[1], coefs1[n]));
// 	mvals[n] = _mm_add_pd(mvals[n], _mm_mul_pd (abc[2], coefs2[n]));
// 	mvals[n] = _mm_add_pd(mvals[n], _mm_mul_pd (abc[3], coefs3[n]));
//       }
// #else
//       for (int k=0; k<4; k++) {
//         __m128d abc = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
//         __m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	
//         for (int n=0; n<Nstop; n+=2) {
//           mvals[n+0] = _mm_add_pd (mvals[n+0], _mm_mul_pd (abc, coefs[n+0]));
//           mvals[n+1] = _mm_add_pd (mvals[n+1], _mm_mul_pd (abc, coefs[n+1]));
// 	}
// 	if (N&1)
// 	  mvals[N-1] = _mm_add_pd (mvals[N-1], _mm_mul_pd (abc, coefs[N-1]));
//       }
// #endif  
//     }
  
//   for (int n=0; n<N; n++)
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
  
// }



// void
// eval_multi_NUBspline_3d_z_vg (multi_NUBspline_3d_z *spline,
// 			     double x, double y, double z,
// 			     complex_double* restrict vals,
// 			     complex_double* restrict grads)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   z -= spline->z_grid.start;
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   double uz = z*spline->z_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
//   double ipartx, iparty, ipartz, tx, ty, tz;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
//   tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int zs = spline->z_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
//     a01  ,   b01,   c01,   a23,    b23,   c23,  
//     da01 ,  db01,  dc01,  da23,   db23,  dc23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);
//   tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
//   tpz23 = _mm_set_pd (tz, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);


//   // z-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);


//   // Zero-out values
//   __m128d mvals[N], mgrads[3*N];
//   for (int n=0; n<N; n++) {
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
//     for (int i=0; i<3; i++) 
//       mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
//   }   

//   __m128d a[4], b[4], c[4], da[4], db[4], dc[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01); 
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01); 
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23); 
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23); 
				 				   
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01); 
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01); 
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23); 
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23); 
				 				   
//   c[0]=_mm_unpacklo_pd(c01,c01); dc[0]=_mm_unpacklo_pd(dc01,dc01); 
//   c[1]=_mm_unpackhi_pd(c01,c01); dc[1]=_mm_unpackhi_pd(dc01,dc01); 
//   c[2]=_mm_unpacklo_pd(c23,c23); dc[2]=_mm_unpacklo_pd(dc23,dc23); 
//   c[3]=_mm_unpackhi_pd(c23,c23); dc[3]=_mm_unpackhi_pd(dc23,dc23); 

//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) 
//       for (int k=0; k<4; k++) {
// 	__m128d abc, d_abc[3];
	
// 	abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
// 	d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
// 	d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
// 	d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);

// 	__m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);

// 	for (int n=0; n<N; n++) {
// 	  mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
// 	  mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
// 	  mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
// 	  mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
// 	}
//       }
  
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   double dzInv = spline->z_grid.delta_inv; 
  
//   for (int n=0; n<N; n++) {
//     complex_double lapl3[3];
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
//     _mm_storeu_pd((double*)(grads+3*n+0), mgrads[3*n+0]);
//     _mm_storeu_pd((double*)(grads+3*n+1), mgrads[3*n+1]);
//     _mm_storeu_pd((double*)(grads+3*n+2), mgrads[3*n+2]);
   
//     grads[3*n+0] *= dxInv;
//     grads[3*n+1] *= dyInv;
//     grads[3*n+2] *= dzInv;
//   }
// }



// void
// eval_multi_NUBspline_3d_z_vgl (multi_NUBspline_3d_z *spline,
// 			      double x, double y, double z,
// 			      complex_double* restrict vals,
// 			      complex_double* restrict grads,
// 			      complex_double* restrict lapl)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   z -= spline->z_grid.start;
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   double uz = z*spline->z_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
//   double ipartx, iparty, ipartz, tx, ty, tz;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
//   tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int zs = spline->z_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
//     a01  ,   b01,   c01,   a23,    b23,   c23,  
//     da01 ,  db01,  dc01,  da23,   db23,  dc23,  
//     d2a01, d2b01, d2c01, d2a23,  d2b23, d2c23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);
//   tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
//   tpz23 = _mm_set_pd (tz, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);


//   // z-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpz01, tpz23, tpz01, tpz23, d2c01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpz01, tpz23, tpz01, tpz23, d2c23);


//   // Zero-out values
//   __m128d mvals[N], mgrads[3*N], mlapl[3*N];
//   for (int n=0; n<N; n++) {
//     mvals[n] = _mm_sub_pd (mvals[n], mvals[n]);
//     for (int i=0; i<3; i++) {
//       mgrads[3*n+i] = _mm_sub_pd (mgrads[3*n+i],mgrads[3*n+i]);
//       mlapl [3*n+i] = _mm_sub_pd (mlapl [3*n+i],mlapl [3*n+i]);
//     }
//   }   

//   __m128d a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01); d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01); d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23); d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23); d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
				 				   				  
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01); d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01); d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23); d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23); d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
				 				   				  
//   c[0]=_mm_unpacklo_pd(c01,c01); dc[0]=_mm_unpacklo_pd(dc01,dc01); d2c[0]=_mm_unpacklo_pd(d2c01,d2c01);
//   c[1]=_mm_unpackhi_pd(c01,c01); dc[1]=_mm_unpackhi_pd(dc01,dc01); d2c[1]=_mm_unpackhi_pd(d2c01,d2c01);
//   c[2]=_mm_unpacklo_pd(c23,c23); dc[2]=_mm_unpacklo_pd(dc23,dc23); d2c[2]=_mm_unpacklo_pd(d2c23,d2c23);
//   c[3]=_mm_unpackhi_pd(c23,c23); dc[3]=_mm_unpackhi_pd(dc23,dc23); d2c[3]=_mm_unpackhi_pd(d2c23,d2c23);

//   // Main computation loop
//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) 
//       for (int k=0; k<4; k++) {
// 	__m128d abc, d_abc[3], d2_abc[3];
	
// 	abc         = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);

// 	d_abc[0]    = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
// 	d_abc[1]    = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
// 	d_abc[2]    = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);

// 	d2_abc[0]   = _mm_mul_pd (_mm_mul_pd(d2a[i],   b[j]),   c[k]);
// 	d2_abc[1]   = _mm_mul_pd (_mm_mul_pd(  a[i], d2b[j]),   c[k]);
// 	d2_abc[2]   = _mm_mul_pd (_mm_mul_pd(  a[i],   b[j]), d2c[k]);
				  

// 	__m128d* restrict coefs = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);

// 	for (int n=0; n<N; n++) {
// 	  mvals[n]      = _mm_add_pd (mvals[n],      _mm_mul_pd (   abc   , coefs[n]));
// 	  mgrads[3*n+0] = _mm_add_pd (mgrads[3*n+0], _mm_mul_pd ( d_abc[0], coefs[n]));
// 	  mgrads[3*n+1] = _mm_add_pd (mgrads[3*n+1], _mm_mul_pd ( d_abc[1], coefs[n]));
// 	  mgrads[3*n+2] = _mm_add_pd (mgrads[3*n+2], _mm_mul_pd ( d_abc[2], coefs[n]));
// 	  mlapl[3*n+0]  = _mm_add_pd (mlapl[3*n+0],  _mm_mul_pd (d2_abc[0], coefs[n]));
// 	  mlapl[3*n+1]  = _mm_add_pd (mlapl[3*n+1],  _mm_mul_pd (d2_abc[1], coefs[n]));
// 	  mlapl[3*n+2]  = _mm_add_pd (mlapl[3*n+2],  _mm_mul_pd (d2_abc[2], coefs[n]));
// 	}
//       }
  
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   double dzInv = spline->z_grid.delta_inv; 
  
//   for (int n=0; n<N; n++) {
//     complex_double lapl3[3];
//     _mm_storeu_pd((double*)(vals+n),mvals[n]);
//     _mm_storeu_pd((double*)(grads+3*n+0), mgrads[3*n+0]);
//     _mm_storeu_pd((double*)(grads+3*n+1), mgrads[3*n+1]);
//     _mm_storeu_pd((double*)(grads+3*n+2), mgrads[3*n+2]);
//     _mm_storeu_pd((double*)(lapl3+0), mlapl[3*n+0]);
//     _mm_storeu_pd((double*)(lapl3+1), mlapl[3*n+1]);
//     _mm_storeu_pd((double*)(lapl3+2), mlapl[3*n+2]);
   
//     grads[3*n+0] *= dxInv;
//     grads[3*n+1] *= dyInv;
//     grads[3*n+2] *= dzInv;
//     lapl3[0] *= dxInv*dxInv;
//     lapl3[1] *= dyInv*dyInv;
//     lapl3[2] *= dzInv*dzInv;
//     lapl[n] = lapl3[0] + lapl3[1] + lapl3[2];
//   }
// }


// void
// eval_multi_NUBspline_3d_z_vgh (multi_NUBspline_3d_z *spline,
// 			      double x, double y, double z,
// 			      complex_double* restrict vals,
// 			      complex_double* restrict grads,
// 			      complex_double* restrict hess)
// {
//   _mm_prefetch ((const char*) &A_d[ 0],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 1],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 2],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 3],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 4],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 5],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 6],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 7],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[ 8],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[ 9],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[10],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[11],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[12],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[13],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[14],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[15],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[16],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[17],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[18],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[19],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[20],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[21],_MM_HINT_T0);  
//   _mm_prefetch ((const char*) &A_d[22],_MM_HINT_T0); _mm_prefetch ((const char*) &A_d[23],_MM_HINT_T0);  

//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;  
//   z -= spline->z_grid.start;
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   double uz = z*spline->z_grid.delta_inv;
//   ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
//   uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
//   uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
//   double ipartx, iparty, ipartz, tx, ty, tz;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
//   tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
//   int xs = spline->x_stride;
//   int ys = spline->y_stride;
//   int zs = spline->z_stride;
//   int N  = spline->num_splines;

//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]

//   // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
//   // A is 4x4 matrix given by the rows A0, A1, A2, A3
//   __m128d tpx01, tpx23, tpy01, tpy23, tpz01, tpz23,
//     a01  ,   b01,   c01,   a23,    b23,   c23,  
//     da01 ,  db01,  dc01,  da23,   db23,  dc23,  
//     d2a01, d2b01, d2c01, d2a23,  d2b23, d2c23;

//   tpx01 = _mm_set_pd (tx*tx*tx, tx*tx);
//   tpx23 = _mm_set_pd (tx, 1.0);
//   tpy01 = _mm_set_pd (ty*ty*ty, ty*ty);
//   tpy23 = _mm_set_pd (ty, 1.0);
//   tpz01 = _mm_set_pd (tz*tz*tz, tz*tz);
//   tpz23 = _mm_set_pd (tz, 1.0);

//   // x-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpx01, tpx23, tpx01, tpx23,   a01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpx01, tpx23, tpx01, tpx23,   a23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpx01, tpx23, tpx01, tpx23,  da01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpx01, tpx23, tpx01, tpx23,  da23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpx01, tpx23, tpx01, tpx23, d2a01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpx01, tpx23, tpx01, tpx23, d2a23);

//   // y-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpy01, tpy23, tpy01, tpy23,   b01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpy01, tpy23, tpy01, tpy23,   b23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpy01, tpy23, tpy01, tpy23,  db01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpy01, tpy23, tpy01, tpy23,  db23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpy01, tpy23, tpy01, tpy23, d2b01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpy01, tpy23, tpy01, tpy23, d2b23);


//   // z-dependent vectors
//   _MM_DDOT4_PD (A_d[ 0], A_d[ 1], A_d[ 2], A_d[ 3], tpz01, tpz23, tpz01, tpz23,   c01);
//   _MM_DDOT4_PD (A_d[ 4], A_d[ 5], A_d[ 6], A_d[ 7], tpz01, tpz23, tpz01, tpz23,   c23);
//   _MM_DDOT4_PD (A_d[ 8], A_d[ 9], A_d[10], A_d[11], tpz01, tpz23, tpz01, tpz23,  dc01);
//   _MM_DDOT4_PD (A_d[12], A_d[13], A_d[14], A_d[15], tpz01, tpz23, tpz01, tpz23,  dc23);
//   _MM_DDOT4_PD (A_d[16], A_d[17], A_d[18], A_d[19], tpz01, tpz23, tpz01, tpz23, d2c01);
//   _MM_DDOT4_PD (A_d[20], A_d[21], A_d[22], A_d[23], tpz01, tpz23, tpz01, tpz23, d2c23);


//   // Zero-out values
//   //__m128d mvals[N], mgrads[3*N], mhess[6*N];
//   __m128d mpack[10*N];
//   for (int n=0; n<10*N; n++) 
//     mpack[n] = _mm_setzero_pd();

//   __m128d a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
//   a[0]=_mm_unpacklo_pd(a01,a01); da[0]=_mm_unpacklo_pd(da01,da01); d2a[0]=_mm_unpacklo_pd(d2a01,d2a01);
//   a[1]=_mm_unpackhi_pd(a01,a01); da[1]=_mm_unpackhi_pd(da01,da01); d2a[1]=_mm_unpackhi_pd(d2a01,d2a01);
//   a[2]=_mm_unpacklo_pd(a23,a23); da[2]=_mm_unpacklo_pd(da23,da23); d2a[2]=_mm_unpacklo_pd(d2a23,d2a23);
//   a[3]=_mm_unpackhi_pd(a23,a23); da[3]=_mm_unpackhi_pd(da23,da23); d2a[3]=_mm_unpackhi_pd(d2a23,d2a23);
				 				   				  
//   b[0]=_mm_unpacklo_pd(b01,b01); db[0]=_mm_unpacklo_pd(db01,db01); d2b[0]=_mm_unpacklo_pd(d2b01,d2b01);
//   b[1]=_mm_unpackhi_pd(b01,b01); db[1]=_mm_unpackhi_pd(db01,db01); d2b[1]=_mm_unpackhi_pd(d2b01,d2b01);
//   b[2]=_mm_unpacklo_pd(b23,b23); db[2]=_mm_unpacklo_pd(db23,db23); d2b[2]=_mm_unpacklo_pd(d2b23,d2b23);
//   b[3]=_mm_unpackhi_pd(b23,b23); db[3]=_mm_unpackhi_pd(db23,db23); d2b[3]=_mm_unpackhi_pd(d2b23,d2b23);
				 				   				  
//   c[0]=_mm_unpacklo_pd(c01,c01); dc[0]=_mm_unpacklo_pd(dc01,dc01); d2c[0]=_mm_unpacklo_pd(d2c01,d2c01);
//   c[1]=_mm_unpackhi_pd(c01,c01); dc[1]=_mm_unpackhi_pd(dc01,dc01); d2c[1]=_mm_unpackhi_pd(d2c01,d2c01);
//   c[2]=_mm_unpacklo_pd(c23,c23); dc[2]=_mm_unpacklo_pd(dc23,dc23); d2c[2]=_mm_unpacklo_pd(d2c23,d2c23);
//   c[3]=_mm_unpackhi_pd(c23,c23); dc[3]=_mm_unpackhi_pd(dc23,dc23); d2c[3]=_mm_unpackhi_pd(d2c23,d2c23);
 
//   // Main computation loop
//   const int bs = 32;
//   for (int nstart=0; nstart<N; nstart += bs) {
//     for (int i=0; i<4; i++)
//       for (int j=0; j<4; j++) {
// 	  __m128d abc[40];
// 	  __m128d* restrict c0 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+0)*zs);
// 	  __m128d* restrict c1 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+1)*zs);
// 	  __m128d* restrict c2 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+2)*zs);
// 	  __m128d* restrict c3 = (__m128d*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+3)*zs);

// #ifdef USE_PREFETCH_VGH
// 	  int nextIndex = i<<4 + j<<2 + k + 1;
// 	  int iNext = nextIndex >> 4; 
// 	  int jNext = (nextIndex >> 2) & 3;
// 	  int kNext = nextIndex & 3;
// 	  if (nextIndex < 64) {
// 	    __m128d* restrict nextCoefs = (__m128d*)(spline->coefs + (ix+iNext)*xs + (iy +jNext)*ys + (iz+kNext)*zs);
// 	    for (int i=0,n=nstart; (n<N && i<bs); n++,i++)
// 	      _mm_prefetch((const char*) &nextCoefs[n], _MM_HINT_NTA);
// 	  }
// #endif
// 	  _mm_prefetch((const char*) &(c0[0]), _MM_HINT_T0);
// 	  _mm_prefetch((const char*) &(c1[0]), _MM_HINT_T0);
// 	  _mm_prefetch((const char*) &(c2[0]), _MM_HINT_T0);
// 	  _mm_prefetch((const char*) &(c3[0]), _MM_HINT_T0);
// 	  _mm_prefetch((const char*) &(c3[0]), _MM_HINT_T0);
// 	  for (int k=0; k<4; k++) {
// 	    abc[k+4*0]   = _mm_mul_pd (_mm_mul_pd(a[i], b[j]), c[k]);
	    
// 	    abc[k+4*1]   = _mm_mul_pd (_mm_mul_pd(da[i],  b[j]),  c[k]);
// 	    abc[k+4*2]   = _mm_mul_pd (_mm_mul_pd( a[i], db[j]),  c[k]);
// 	    abc[k+4*3]   = _mm_mul_pd (_mm_mul_pd( a[i],  b[j]), dc[k]);
	    
// 	    abc[k+4*4]   = _mm_mul_pd (_mm_mul_pd(d2a[i],   b[j]),   c[k]);
// 	    abc[k+4*5]   = _mm_mul_pd (_mm_mul_pd( da[i],  db[j]),   c[k]);
// 	    abc[k+4*6]   = _mm_mul_pd (_mm_mul_pd( da[i],   b[j]),  dc[k]);
// 	    abc[k+4*7]   = _mm_mul_pd (_mm_mul_pd(  a[i], d2b[j]),   c[k]);
// 	    abc[k+4*8]   = _mm_mul_pd (_mm_mul_pd(  a[i],  db[j]),  dc[k]);
// 	    abc[k+4*9]   = _mm_mul_pd (_mm_mul_pd(  a[i],   b[j]), d2c[k]);
// 	  }
// 	  int end; 
// 	  if (N < nstart+bs)   end = N; else end = nstart+bs;
	  
// 	  for (int n=nstart; n<end; n++) 
// 	    for (int s=0; s<10; s++) {
// 	      __m128d p0 = _mm_mul_pd(abc[4*s+0], c0[n]);
// 	      __m128d p1 = _mm_mul_pd(abc[4*s+1], c1[n]);
// 	      __m128d p2 = _mm_mul_pd(abc[4*s+2], c2[n]);
// 	      __m128d p3 = _mm_mul_pd(abc[4*s+3], c3[n]);
// 	      __m128d sum0 = _mm_add_pd (p0, p1);
// 	      __m128d sum1 = _mm_add_pd (p2, p3);
// 	      __m128d sum2 = _mm_add_pd (sum0, sum1);

// 	      mpack[10*n+s] = _mm_add_pd (mpack[10*n+s], sum2);
// 	      //					  mm_add_pd(_mm_mul_pd (   abc[s], c0[n])), _mm_mul_pd(abc[s], c1[n]))
// 	  //mpack[n+s*N] = _mm_add_pd (mpack[n+s*N], _mm_mul_pd (
// 	  //abc[s], coefs[n]));
// 	    }
	  
//       }
//   }
    
//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   double dzInv = spline->z_grid.delta_inv; 
  
//   for (int n=0; n<N; n++) {
//     _mm_storeu_pd((double*)(vals+n)     , mpack[10*n+0]);
//     _mm_storeu_pd((double*)(grads+3*n+0), mpack[10*n+1]);
//     _mm_storeu_pd((double*)(grads+3*n+1), mpack[10*n+2]);
//     _mm_storeu_pd((double*)(grads+3*n+2), mpack[10*n+3]);
 
//     _mm_storeu_pd((double*)(hess+9*n+0),  mpack[10*n+4]);

//     _mm_storeu_pd((double*)(hess+9*n+1),  mpack[10*n+5]);
//     _mm_storeu_pd((double*)(hess+9*n+3),  mpack[10*n+5]);

//     _mm_storeu_pd((double*)(hess+9*n+2),  mpack[10*n+6]);
//     _mm_storeu_pd((double*)(hess+9*n+6),  mpack[10*n+6]);

//     _mm_storeu_pd((double*)(hess+9*n+4),  mpack[10*n+7]);

//     _mm_storeu_pd((double*)(hess+9*n+5),  mpack[10*n+8]);
//     _mm_storeu_pd((double*)(hess+9*n+7),  mpack[10*n+8]);

//     _mm_storeu_pd((double*)(hess+9*n+8),  mpack[10*n+9]);
//   }
//   for (int n=0; n<N; n++) {
//     grads[3*n+0] *= dxInv;
//     grads[3*n+1] *= dyInv;
//     grads[3*n+2] *= dzInv;
//     hess[9*n+0]  *= dxInv*dxInv;
//     hess[9*n+4]  *= dyInv*dyInv;
//     hess[9*n+8]  *= dzInv*dzInv;

//     hess[9*n+1]  *= dxInv*dyInv;
//     hess[9*n+3]  *= dxInv*dyInv;

//     hess[9*n+2]  *= dxInv*dzInv;
//     hess[9*n+6]  *= dxInv*dzInv;

//     hess[9*n+5]  *= dyInv*dzInv;
//     hess[9*n+7]  *= dyInv*dzInv;
// //     // Copy hessian elements into lower half of 3x3 matrix
// //     hess[9*n+3] = hess[9*n+1];
// //     hess[9*n+6] = hess[9*n+2];
// //     hess[9*n+7] = hess[9*n+5];
//   }
// }

#endif
