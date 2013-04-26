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

#ifndef BSPLINE_EVAL_STD_C_H
#define BSPLINE_EVAL_STD_C_H

#include <math.h>
#include <stdio.h>

extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_c (UBspline_1d_c * restrict spline,
                    double x, float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float tp[4], b[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;
  // Basis functions
  b[0] = (Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3]);
  b[1] = (Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3]);
  b[2] = (Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3]);
  b[3] = (Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]);
#define Pr(i) coefs[2*i]
#define Pi(i) coefs[2*i+1]
  // Real part
  val[0] = Pr(i+0)*b[0] + Pr(i+1)*b[1] + Pr(i+2)*b[2] + Pr(i+3)*b[3];
  // Imaginary part
  val[1] = Pi(i+0)*b[0] + Pi(i+1)*b[1] + Pi(i+2)*b[2] + Pi(i+3)*b[3];
#undef Pr
#undef Pi
}

/* Value and first derivative */
inline void
eval_UBspline_1d_c_vg (UBspline_1d_c * restrict spline, double x,
                       float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float tp[4], b[4], db[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;
  // Basis functions
  b[0] = (Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3]);
  b[1] = (Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3]);
  b[2] = (Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3]);
  b[3] = (Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]);
  // 1st derivative
  db[0] = (dAf[ 0]*tp[0] + dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3]);
  db[1] = (dAf[ 4]*tp[0] + dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3]);
  db[2] = (dAf[ 8]*tp[0] + dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3]);
  db[3] = (dAf[12]*tp[0] + dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]);
  // Real part
  val[0]  = Pr(i+0)* b[0] + Pr(i+1)* b[1] + Pr(i+2)* b[2] + Pr(i+3)* b[3];
  grad[0] = Pr(i+0)*db[0] + Pr(i+1)*db[1] + Pr(i+2)*db[2] + Pr(i+3)*db[3];
  // Imaginary part
  val[1]  = Pi(i+0)* b[0] + Pi(i+1)* b[1] + Pi(i+2)* b[2] + Pi(i+3)* b[3];
  grad[1] = Pi(i+0)*db[0] + Pi(i+1)*db[1] + Pi(i+2)*db[2] + Pi(i+3)*db[3];
}

/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_c_vgl (UBspline_1d_c * restrict spline, double x,
                        float* restrict val, float* restrict grad,
                        float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float tp[4], b[4], db[4], d2b[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;
  // Basis functions
  b[0] = (Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3]);
  b[1] = (Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3]);
  b[2] = (Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3]);
  b[3] = (Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]);
  // 1st derivative
  db[0] = (dAf[ 0]*tp[0] + dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3]);
  db[1] = (dAf[ 4]*tp[0] + dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3]);
  db[2] = (dAf[ 8]*tp[0] + dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3]);
  db[3] = (dAf[12]*tp[0] + dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]);
  // 2nd derivative
  d2b[0] = (d2Af[ 0]*tp[0] + d2Af[ 1]*tp[1] + d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3]);
  d2b[1] = (d2Af[ 4]*tp[0] + d2Af[ 5]*tp[1] + d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3]);
  d2b[2] = (d2Af[ 8]*tp[0] + d2Af[ 9]*tp[1] + d2Af[10]*tp[2] + d2Af[11]*tp[3]);
  d2b[3] = (d2Af[12]*tp[0] + d2Af[13]*tp[1] + d2Af[14]*tp[2] + d2Af[15]*tp[3]);
  // Real part
  val[0]  = Pr(i+0)*  b[0] + Pr(i+1)*  b[1] + Pr(i+2)*  b[2] + Pr(i+3)*  b[3];
  grad[0] = Pr(i+0)* db[0] + Pr(i+1)* db[1] + Pr(i+2)* db[2] + Pr(i+3)* db[3];
  lapl[0] = Pr(i+0)*d2b[0] + Pr(i+1)*d2b[1] + Pr(i+2)*d2b[2] + Pr(i+3)*d2b[3];
  // Imaginary part
  val[1]  = Pi(i+0)*  b[0] + Pi(i+1)*  b[1] + Pi(i+2)*  b[2] + Pi(i+3)*  b[3];
  grad[1] = Pi(i+0)* db[0] + Pi(i+1)* db[1] + Pi(i+2)* db[2] + Pi(i+3)* db[3];
  lapl[1] = Pi(i+0)*d2b[0] + Pi(i+1)*d2b[1] + Pi(i+2)*d2b[2] + Pi(i+3)*d2b[3];
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_c (UBspline_2d_c * restrict spline,
                    double x, double y, float* restrict val)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  float tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  int xs = spline->x_ctride;
#define Pr(i,j) coefs[2*((ix+(i))*xs+iy+(j))]
#define Pi(i,j) coefs[2*((ix+(i))*xs+iy+(j))+1]
  val[0] = (a[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
            a[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
            a[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
            a[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  val[1] = (a[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
            a[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
            a[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
            a[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
#undef Pr
#undef Pi
}


/* Value and gradient */
inline void
eval_UBspline_2d_c_vg (UBspline_2d_c * restrict spline,
                       double x, double y,
                       float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  b[0]  = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  int xs = spline->x_ctride;
#define Pr(i,j) coefs[2*((ix+(i))*xs+iy+(j))]
#define Pi(i,j) coefs[2*((ix+(i))*xs+iy+(j))+1]
  // Real part
  val[0] =
    (a[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
     a[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
     a[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
     a[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  val[1] =
    (a[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
     a[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
     a[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
     a[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[0] = spline->x_grid.delta_inv *
            (da[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
             da[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
             da[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
             da[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  grad[1] = spline->x_grid.delta_inv *
            (da[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
             da[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
             da[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
             da[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[2] = spline->y_grid.delta_inv *
            (a[0]*(Pr(0,0)*db[0]+Pr(0,1)*db[1]+Pr(0,2)*db[2]+Pr(0,3)*db[3])+
             a[1]*(Pr(1,0)*db[0]+Pr(1,1)*db[1]+Pr(1,2)*db[2]+Pr(1,3)*db[3])+
             a[2]*(Pr(2,0)*db[0]+Pr(2,1)*db[1]+Pr(2,2)*db[2]+Pr(2,3)*db[3])+
             a[3]*(Pr(3,0)*db[0]+Pr(3,1)*db[1]+Pr(3,2)*db[2]+Pr(3,3)*db[3]));
  // Imag part
  grad[3] = spline->y_grid.delta_inv *
            (a[0]*(Pi(0,0)*db[0]+Pi(0,1)*db[1]+Pi(0,2)*db[2]+Pi(0,3)*db[3])+
             a[1]*(Pi(1,0)*db[0]+Pi(1,1)*db[1]+Pi(1,2)*db[2]+Pi(1,3)*db[3])+
             a[2]*(Pi(2,0)*db[0]+Pi(2,1)*db[1]+Pi(2,2)*db[2]+Pi(2,3)*db[3])+
             a[3]*(Pi(3,0)*db[0]+Pi(3,1)*db[1]+Pi(3,2)*db[2]+Pi(3,3)*db[3]));
#undef Pr
#undef Pi
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_c_vgl (UBspline_2d_c * restrict spline,
                        double x, double y, float* restrict val,
                        float* restrict grad, float* restrict lapl)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);
  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  int xs = spline->x_ctride;
#define Pr(i,j) coefs[2*((ix+(i))*xs+iy+(j))]
#define Pi(i,j) coefs[2*((ix+(i))*xs+iy+(j))+1]
  // Real part
  val[0] =
    (a[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
     a[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
     a[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
     a[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  val[1] =
    (a[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
     a[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
     a[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
     a[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[0] = spline->x_grid.delta_inv *
            (da[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
             da[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
             da[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
             da[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  grad[1] = spline->x_grid.delta_inv *
            (da[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
             da[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
             da[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
             da[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[2] = spline->y_grid.delta_inv *
            (a[0]*(Pr(0,0)*db[0]+Pr(0,1)*db[1]+Pr(0,2)*db[2]+Pr(0,3)*db[3])+
             a[1]*(Pr(1,0)*db[0]+Pr(1,1)*db[1]+Pr(1,2)*db[2]+Pr(1,3)*db[3])+
             a[2]*(Pr(2,0)*db[0]+Pr(2,1)*db[1]+Pr(2,2)*db[2]+Pr(2,3)*db[3])+
             a[3]*(Pr(3,0)*db[0]+Pr(3,1)*db[1]+Pr(3,2)*db[2]+Pr(3,3)*db[3]));
  // Imag part
  grad[3] = spline->y_grid.delta_inv *
            (a[0]*(Pi(0,0)*db[0]+Pi(0,1)*db[1]+Pi(0,2)*db[2]+Pi(0,3)*db[3])+
             a[1]*(Pi(1,0)*db[0]+Pi(1,1)*db[1]+Pi(1,2)*db[2]+Pi(1,3)*db[3])+
             a[2]*(Pi(2,0)*db[0]+Pi(2,1)*db[1]+Pi(2,2)*db[2]+Pi(2,3)*db[3])+
             a[3]*(Pi(3,0)*db[0]+Pi(3,1)*db[1]+Pi(3,2)*db[2]+Pi(3,3)*db[3]));
  // Real Part
  lapl[0]   =
    spline->x_grid.delta_inv * spline->x_grid.delta_inv *
    (a[0]*(Pr(0,0)*d2b[0]+Pr(0,1)*d2b[1]+Pr(0,2)*d2b[2]+Pr(0,3)*d2b[3])+
     a[1]*(Pr(1,0)*d2b[0]+Pr(1,1)*d2b[1]+Pr(1,2)*d2b[2]+Pr(1,3)*d2b[3])+
     a[2]*(Pr(2,0)*d2b[0]+Pr(2,1)*d2b[1]+Pr(2,2)*d2b[2]+Pr(2,3)*d2b[3])+
     a[3]*(Pr(3,0)*d2b[0]+Pr(3,1)*d2b[1]+Pr(3,2)*d2b[2]+Pr(3,3)*d2b[3])) +
    spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (d2a[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
     d2a[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
     d2a[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
     d2a[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  lapl[1]   =
    spline->x_grid.delta_inv * spline->x_grid.delta_inv *
    (a[0]*(Pi(0,0)*d2b[0]+Pi(0,1)*d2b[1]+Pi(0,2)*d2b[2]+Pi(0,3)*d2b[3])+
     a[1]*(Pi(1,0)*d2b[0]+Pi(1,1)*d2b[1]+Pi(1,2)*d2b[2]+Pi(1,3)*d2b[3])+
     a[2]*(Pi(2,0)*d2b[0]+Pi(2,1)*d2b[1]+Pi(2,2)*d2b[2]+Pi(2,3)*d2b[3])+
     a[3]*(Pi(3,0)*d2b[0]+Pi(3,1)*d2b[1]+Pi(3,2)*d2b[2]+Pi(3,3)*d2b[3])) +
    spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (d2a[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
     d2a[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
     d2a[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
     d2a[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
#undef Pr
#undef Pi
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_c_vgh (UBspline_2d_c * restrict spline,
                        double x, double y, float* restrict val,
                        float* restrict grad, float* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);
  b[0]   = (  Af[ 0]*tpy[0] +   Af[ 1]*tpy[1] +  Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]   = (  Af[ 4]*tpy[0] +   Af[ 5]*tpy[1] +  Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]   = (  Af[ 8]*tpy[0] +   Af[ 9]*tpy[1] +  Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]   = (  Af[12]*tpy[0] +   Af[13]*tpy[1] +  Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0]  = ( dAf[ 1]*tpy[1] +  dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1]  = ( dAf[ 5]*tpy[1] +  dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2]  = ( dAf[ 9]*tpy[1] +  dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3]  = ( dAf[13]*tpy[1] +  dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  int xs = spline->x_ctride;
#define Pr(i,j) coefs[2*((ix+(i))*xs+iy+(j))]
#define Pi(i,j) coefs[2*((ix+(i))*xs+iy+(j))+1]
  // Real part
  val[0] =
    (a[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
     a[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
     a[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
     a[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  val[1] =
    (a[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
     a[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
     a[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
     a[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[0] = spline->x_grid.delta_inv *
            (da[0]*(Pr(0,0)*b[0]+Pr(0,1)*b[1]+Pr(0,2)*b[2]+Pr(0,3)*b[3])+
             da[1]*(Pr(1,0)*b[0]+Pr(1,1)*b[1]+Pr(1,2)*b[2]+Pr(1,3)*b[3])+
             da[2]*(Pr(2,0)*b[0]+Pr(2,1)*b[1]+Pr(2,2)*b[2]+Pr(2,3)*b[3])+
             da[3]*(Pr(3,0)*b[0]+Pr(3,1)*b[1]+Pr(3,2)*b[2]+Pr(3,3)*b[3]));
  // Imag part
  grad[1] = spline->x_grid.delta_inv *
            (da[0]*(Pi(0,0)*b[0]+Pi(0,1)*b[1]+Pi(0,2)*b[2]+Pi(0,3)*b[3])+
             da[1]*(Pi(1,0)*b[0]+Pi(1,1)*b[1]+Pi(1,2)*b[2]+Pi(1,3)*b[3])+
             da[2]*(Pi(2,0)*b[0]+Pi(2,1)*b[1]+Pi(2,2)*b[2]+Pi(2,3)*b[3])+
             da[3]*(Pi(3,0)*b[0]+Pi(3,1)*b[1]+Pi(3,2)*b[2]+Pi(3,3)*b[3]));
  // Real part
  grad[2] = spline->y_grid.delta_inv *
            (a[0]*(Pr(0,0)*db[0]+Pr(0,1)*db[1]+Pr(0,2)*db[2]+Pr(0,3)*db[3])+
             a[1]*(Pr(1,0)*db[0]+Pr(1,1)*db[1]+Pr(1,2)*db[2]+Pr(1,3)*db[3])+
             a[2]*(Pr(2,0)*db[0]+Pr(2,1)*db[1]+Pr(2,2)*db[2]+Pr(2,3)*db[3])+
             a[3]*(Pr(3,0)*db[0]+Pr(3,1)*db[1]+Pr(3,2)*db[2]+Pr(3,3)*db[3]));
  // Imag part
  grad[3] = spline->y_grid.delta_inv *
            (a[0]*(Pi(0,0)*db[0]+Pi(0,1)*db[1]+Pi(0,2)*db[2]+Pi(0,3)*db[3])+
             a[1]*(Pi(1,0)*db[0]+Pi(1,1)*db[1]+Pi(1,2)*db[2]+Pi(1,3)*db[3])+
             a[2]*(Pi(2,0)*db[0]+Pi(2,1)*db[1]+Pi(2,2)*db[2]+Pi(2,3)*db[3])+
             a[3]*(Pi(3,0)*db[0]+Pi(3,1)*db[1]+Pi(3,2)*db[2]+Pi(3,3)*db[3]));
  // Real part
  hess[0] = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
            (d2a[0]*(Pr(0,0)*  b[0]+Pr(0,1)*  b[1]+Pr(0,2)*  b[2]+Pr(0,3)*  b[3])+
             d2a[1]*(Pr(1,0)*  b[0]+Pr(1,1)*  b[1]+Pr(1,2)*  b[2]+Pr(1,3)*  b[3])+
             d2a[2]*(Pr(2,0)*  b[0]+Pr(2,1)*  b[1]+Pr(2,2)*  b[2]+Pr(2,3)*  b[3])+
             d2a[3]*(Pr(3,0)*  b[0]+Pr(3,1)*  b[1]+Pr(3,2)*  b[2]+Pr(3,3)*  b[3]));
  // Imag part
  hess[1] = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
            (d2a[0]*(Pi(0,0)*  b[0]+Pi(0,1)*  b[1]+Pi(0,2)*  b[2]+Pi(0,3)*  b[3])+
             d2a[1]*(Pi(1,0)*  b[0]+Pi(1,1)*  b[1]+Pi(1,2)*  b[2]+Pi(1,3)*  b[3])+
             d2a[2]*(Pi(2,0)*  b[0]+Pi(2,1)*  b[1]+Pi(2,2)*  b[2]+Pi(2,3)*  b[3])+
             d2a[3]*(Pi(3,0)*  b[0]+Pi(3,1)*  b[1]+Pi(3,2)*  b[2]+Pi(3,3)*  b[3]));
  // Real part
  hess[2] = spline->x_grid.delta_inv * spline->y_grid.delta_inv *
            ( da[0]*(Pr(0,0)* db[0]+Pr(0,1)* db[1]+Pr(0,2)* db[2]+Pr(0,3)* db[3])+
              da[1]*(Pr(1,0)* db[0]+Pr(1,1)* db[1]+Pr(1,2)* db[2]+Pr(1,3)* db[3])+
              da[2]*(Pr(2,0)* db[0]+Pr(2,1)* db[1]+Pr(2,2)* db[2]+Pr(2,3)* db[3])+
              da[3]*(Pr(3,0)* db[0]+Pr(3,1)* db[1]+Pr(3,2)* db[2]+Pr(3,3)* db[3]));
  // Imag part
  hess[3] = spline->x_grid.delta_inv * spline->y_grid.delta_inv *
            ( da[0]*(Pi(0,0)* db[0]+Pi(0,1)* db[1]+Pi(0,2)* db[2]+Pi(0,3)* db[3])+
              da[1]*(Pi(1,0)* db[0]+Pi(1,1)* db[1]+Pi(1,2)* db[2]+Pi(1,3)* db[3])+
              da[2]*(Pi(2,0)* db[0]+Pi(2,1)* db[1]+Pi(2,2)* db[2]+Pi(2,3)* db[3])+
              da[3]*(Pi(3,0)* db[0]+Pi(3,1)* db[1]+Pi(3,2)* db[2]+Pi(3,3)* db[3]));
  // Real part
  hess[6] = spline->y_grid.delta_inv * spline->y_grid.delta_inv *
            (  a[0]*(Pr(0,0)*d2b[0]+Pr(0,1)*d2b[1]+Pr(0,2)*d2b[2]+Pr(0,3)*d2b[3])+
               a[1]*(Pr(1,0)*d2b[0]+Pr(1,1)*d2b[1]+Pr(1,2)*d2b[2]+Pr(1,3)*d2b[3])+
               a[2]*(Pr(2,0)*d2b[0]+Pr(2,1)*d2b[1]+Pr(2,2)*d2b[2]+Pr(2,3)*d2b[3])+
               a[3]*(Pr(3,0)*d2b[0]+Pr(3,1)*d2b[1]+Pr(3,2)*d2b[2]+Pr(3,3)*d2b[3]));
  // Imag part
  hess[7] = spline->y_grid.delta_inv * spline->y_grid.delta_inv *
            (  a[0]*(Pi(0,0)*d2b[0]+Pi(0,1)*d2b[1]+Pi(0,2)*d2b[2]+Pi(0,3)*d2b[3])+
               a[1]*(Pi(1,0)*d2b[0]+Pi(1,1)*d2b[1]+Pi(1,2)*d2b[2]+Pi(1,3)*d2b[3])+
               a[2]*(Pi(2,0)*d2b[0]+Pi(2,1)*d2b[1]+Pi(2,2)*d2b[2]+Pi(2,3)*d2b[3])+
               a[3]*(Pi(3,0)*d2b[0]+Pi(3,1)*d2b[1]+Pi(3,2)*d2b[2]+Pi(3,3)*d2b[3]));
  // Real part
  hess[4] = hess[2];
  // Imag part
  hess[5] = hess[3];
#undef Pr
#undef Pi
}


/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_c (UBspline_3d_c * restrict spline,
                    double x, double y, double z,
                    float* restrict val)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modff (uy, &iparty);
  int iy = (int) iparty;
  tz = modff (uz, &ipartz);
  int iz = (int) ipartz;
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;
  tpz[1] = tz*tz;
  tpz[2] = tz;
  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  int xs = spline->x_ctride;
  int ys = spline->y_ctride;
#define Pr(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))]
#define Pi(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))+1]
  val[0] = (a[0]*(b[0]*(Pr(0,0,0)*c[0]+Pr(0,0,1)*c[1]+Pr(0,0,2)*c[2]+Pr(0,0,3)*c[3])+
                  b[1]*(Pr(0,1,0)*c[0]+Pr(0,1,1)*c[1]+Pr(0,1,2)*c[2]+Pr(0,1,3)*c[3])+
                  b[2]*(Pr(0,2,0)*c[0]+Pr(0,2,1)*c[1]+Pr(0,2,2)*c[2]+Pr(0,2,3)*c[3])+
                  b[3]*(Pr(0,3,0)*c[0]+Pr(0,3,1)*c[1]+Pr(0,3,2)*c[2]+Pr(0,3,3)*c[3]))+
            a[1]*(b[0]*(Pr(1,0,0)*c[0]+Pr(1,0,1)*c[1]+Pr(1,0,2)*c[2]+Pr(1,0,3)*c[3])+
                  b[1]*(Pr(1,1,0)*c[0]+Pr(1,1,1)*c[1]+Pr(1,1,2)*c[2]+Pr(1,1,3)*c[3])+
                  b[2]*(Pr(1,2,0)*c[0]+Pr(1,2,1)*c[1]+Pr(1,2,2)*c[2]+Pr(1,2,3)*c[3])+
                  b[3]*(Pr(1,3,0)*c[0]+Pr(1,3,1)*c[1]+Pr(1,3,2)*c[2]+Pr(1,3,3)*c[3]))+
            a[2]*(b[0]*(Pr(2,0,0)*c[0]+Pr(2,0,1)*c[1]+Pr(2,0,2)*c[2]+Pr(2,0,3)*c[3])+
                  b[1]*(Pr(2,1,0)*c[0]+Pr(2,1,1)*c[1]+Pr(2,1,2)*c[2]+Pr(2,1,3)*c[3])+
                  b[2]*(Pr(2,2,0)*c[0]+Pr(2,2,1)*c[1]+Pr(2,2,2)*c[2]+Pr(2,2,3)*c[3])+
                  b[3]*(Pr(2,3,0)*c[0]+Pr(2,3,1)*c[1]+Pr(2,3,2)*c[2]+Pr(2,3,3)*c[3]))+
            a[3]*(b[0]*(Pr(3,0,0)*c[0]+Pr(3,0,1)*c[1]+Pr(3,0,2)*c[2]+Pr(3,0,3)*c[3])+
                  b[1]*(Pr(3,1,0)*c[0]+Pr(3,1,1)*c[1]+Pr(3,1,2)*c[2]+Pr(3,1,3)*c[3])+
                  b[2]*(Pr(3,2,0)*c[0]+Pr(3,2,1)*c[1]+Pr(3,2,2)*c[2]+Pr(3,2,3)*c[3])+
                  b[3]*(Pr(3,3,0)*c[0]+Pr(3,3,1)*c[1]+Pr(3,3,2)*c[2]+Pr(3,3,3)*c[3])));
  val[1] = (a[0]*(b[0]*(Pi(0,0,0)*c[0]+Pi(0,0,1)*c[1]+Pi(0,0,2)*c[2]+Pi(0,0,3)*c[3])+
                  b[1]*(Pi(0,1,0)*c[0]+Pi(0,1,1)*c[1]+Pi(0,1,2)*c[2]+Pi(0,1,3)*c[3])+
                  b[2]*(Pi(0,2,0)*c[0]+Pi(0,2,1)*c[1]+Pi(0,2,2)*c[2]+Pi(0,2,3)*c[3])+
                  b[3]*(Pi(0,3,0)*c[0]+Pi(0,3,1)*c[1]+Pi(0,3,2)*c[2]+Pi(0,3,3)*c[3]))+
            a[1]*(b[0]*(Pi(1,0,0)*c[0]+Pi(1,0,1)*c[1]+Pi(1,0,2)*c[2]+Pi(1,0,3)*c[3])+
                  b[1]*(Pi(1,1,0)*c[0]+Pi(1,1,1)*c[1]+Pi(1,1,2)*c[2]+Pi(1,1,3)*c[3])+
                  b[2]*(Pi(1,2,0)*c[0]+Pi(1,2,1)*c[1]+Pi(1,2,2)*c[2]+Pi(1,2,3)*c[3])+
                  b[3]*(Pi(1,3,0)*c[0]+Pi(1,3,1)*c[1]+Pi(1,3,2)*c[2]+Pi(1,3,3)*c[3]))+
            a[2]*(b[0]*(Pi(2,0,0)*c[0]+Pi(2,0,1)*c[1]+Pi(2,0,2)*c[2]+Pi(2,0,3)*c[3])+
                  b[1]*(Pi(2,1,0)*c[0]+Pi(2,1,1)*c[1]+Pi(2,1,2)*c[2]+Pi(2,1,3)*c[3])+
                  b[2]*(Pi(2,2,0)*c[0]+Pi(2,2,1)*c[1]+Pi(2,2,2)*c[2]+Pi(2,2,3)*c[3])+
                  b[3]*(Pi(2,3,0)*c[0]+Pi(2,3,1)*c[1]+Pi(2,3,2)*c[2]+Pi(2,3,3)*c[3]))+
            a[3]*(b[0]*(Pi(3,0,0)*c[0]+Pi(3,0,1)*c[1]+Pi(3,0,2)*c[2]+Pi(3,0,3)*c[3])+
                  b[1]*(Pi(3,1,0)*c[0]+Pi(3,1,1)*c[1]+Pi(3,1,2)*c[2]+Pi(3,1,3)*c[3])+
                  b[2]*(Pi(3,2,0)*c[0]+Pi(3,2,1)*c[1]+Pi(3,2,2)*c[2]+Pi(3,2,3)*c[3])+
                  b[3]*(Pi(3,3,0)*c[0]+Pi(3,3,1)*c[1]+Pi(3,3,2)*c[2]+Pi(3,3,3)*c[3])));
#undef Pr
#undef Pi
}

/* Value and gradient */
inline void
eval_UBspline_3d_c_vg (UBspline_3d_c * restrict spline,
                       double x, double y, double z,
                       float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modff (uy, &iparty);
  int iy = (int) iparty;
  tz = modff (uz, &ipartz);
  int iz = (int) ipartz;
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4],
        cP[32], dcP[32], bcP[8], dbcP[8];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;
  tpz[1] = tz*tz;
  tpz[2] = tz;
  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  int xs = spline->x_ctride;
  int ys = spline->y_ctride;
#define Pr(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))]
#define Pi(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))+1]
  cP[ 0] = (Pr(0,0,0)*c[0]+Pr(0,0,1)*c[1]+Pr(0,0,2)*c[2]+Pr(0,0,3)*c[3]);
  cP[ 1] = (Pi(0,0,0)*c[0]+Pi(0,0,1)*c[1]+Pi(0,0,2)*c[2]+Pi(0,0,3)*c[3]);
  cP[ 2] = (Pr(0,1,0)*c[0]+Pr(0,1,1)*c[1]+Pr(0,1,2)*c[2]+Pr(0,1,3)*c[3]);
  cP[ 3] = (Pi(0,1,0)*c[0]+Pi(0,1,1)*c[1]+Pi(0,1,2)*c[2]+Pi(0,1,3)*c[3]);
  cP[ 4] = (Pr(0,2,0)*c[0]+Pr(0,2,1)*c[1]+Pr(0,2,2)*c[2]+Pr(0,2,3)*c[3]);
  cP[ 5] = (Pi(0,2,0)*c[0]+Pi(0,2,1)*c[1]+Pi(0,2,2)*c[2]+Pi(0,2,3)*c[3]);
  cP[ 6] = (Pr(0,3,0)*c[0]+Pr(0,3,1)*c[1]+Pr(0,3,2)*c[2]+Pr(0,3,3)*c[3]);
  cP[ 7] = (Pi(0,3,0)*c[0]+Pi(0,3,1)*c[1]+Pi(0,3,2)*c[2]+Pi(0,3,3)*c[3]);
  cP[ 8] = (Pr(1,0,0)*c[0]+Pr(1,0,1)*c[1]+Pr(1,0,2)*c[2]+Pr(1,0,3)*c[3]);
  cP[ 9] = (Pi(1,0,0)*c[0]+Pi(1,0,1)*c[1]+Pi(1,0,2)*c[2]+Pi(1,0,3)*c[3]);
  cP[10] = (Pr(1,1,0)*c[0]+Pr(1,1,1)*c[1]+Pr(1,1,2)*c[2]+Pr(1,1,3)*c[3]);
  cP[11] = (Pi(1,1,0)*c[0]+Pi(1,1,1)*c[1]+Pi(1,1,2)*c[2]+Pi(1,1,3)*c[3]);
  cP[12] = (Pr(1,2,0)*c[0]+Pr(1,2,1)*c[1]+Pr(1,2,2)*c[2]+Pr(1,2,3)*c[3]);
  cP[13] = (Pi(1,2,0)*c[0]+Pi(1,2,1)*c[1]+Pi(1,2,2)*c[2]+Pi(1,2,3)*c[3]);
  cP[14] = (Pr(1,3,0)*c[0]+Pr(1,3,1)*c[1]+Pr(1,3,2)*c[2]+Pr(1,3,3)*c[3]);
  cP[15] = (Pi(1,3,0)*c[0]+Pi(1,3,1)*c[1]+Pi(1,3,2)*c[2]+Pi(1,3,3)*c[3]);
  cP[16] = (Pr(2,0,0)*c[0]+Pr(2,0,1)*c[1]+Pr(2,0,2)*c[2]+Pr(2,0,3)*c[3]);
  cP[17] = (Pi(2,0,0)*c[0]+Pi(2,0,1)*c[1]+Pi(2,0,2)*c[2]+Pi(2,0,3)*c[3]);
  cP[18] = (Pr(2,1,0)*c[0]+Pr(2,1,1)*c[1]+Pr(2,1,2)*c[2]+Pr(2,1,3)*c[3]);
  cP[19] = (Pi(2,1,0)*c[0]+Pi(2,1,1)*c[1]+Pi(2,1,2)*c[2]+Pi(2,1,3)*c[3]);
  cP[20] = (Pr(2,2,0)*c[0]+Pr(2,2,1)*c[1]+Pr(2,2,2)*c[2]+Pr(2,2,3)*c[3]);
  cP[21] = (Pi(2,2,0)*c[0]+Pi(2,2,1)*c[1]+Pi(2,2,2)*c[2]+Pi(2,2,3)*c[3]);
  cP[22] = (Pr(2,3,0)*c[0]+Pr(2,3,1)*c[1]+Pr(2,3,2)*c[2]+Pr(2,3,3)*c[3]);
  cP[23] = (Pi(2,3,0)*c[0]+Pi(2,3,1)*c[1]+Pi(2,3,2)*c[2]+Pi(2,3,3)*c[3]);
  cP[24] = (Pr(3,0,0)*c[0]+Pr(3,0,1)*c[1]+Pr(3,0,2)*c[2]+Pr(3,0,3)*c[3]);
  cP[25] = (Pi(3,0,0)*c[0]+Pi(3,0,1)*c[1]+Pi(3,0,2)*c[2]+Pi(3,0,3)*c[3]);
  cP[26] = (Pr(3,1,0)*c[0]+Pr(3,1,1)*c[1]+Pr(3,1,2)*c[2]+Pr(3,1,3)*c[3]);
  cP[27] = (Pi(3,1,0)*c[0]+Pi(3,1,1)*c[1]+Pi(3,1,2)*c[2]+Pi(3,1,3)*c[3]);
  cP[28] = (Pr(3,2,0)*c[0]+Pr(3,2,1)*c[1]+Pr(3,2,2)*c[2]+Pr(3,2,3)*c[3]);
  cP[29] = (Pi(3,2,0)*c[0]+Pi(3,2,1)*c[1]+Pi(3,2,2)*c[2]+Pi(3,2,3)*c[3]);
  cP[30] = (Pr(3,3,0)*c[0]+Pr(3,3,1)*c[1]+Pr(3,3,2)*c[2]+Pr(3,3,3)*c[3]);
  cP[31] = (Pi(3,3,0)*c[0]+Pi(3,3,1)*c[1]+Pi(3,3,2)*c[2]+Pi(3,3,3)*c[3]);
  dcP[ 0] = (Pr(0,0,0)*dc[0]+Pr(0,0,1)*dc[1]+Pr(0,0,2)*dc[2]+Pr(0,0,3)*dc[3]);
  dcP[ 1] = (Pi(0,0,0)*dc[0]+Pi(0,0,1)*dc[1]+Pi(0,0,2)*dc[2]+Pi(0,0,3)*dc[3]);
  dcP[ 2] = (Pr(0,1,0)*dc[0]+Pr(0,1,1)*dc[1]+Pr(0,1,2)*dc[2]+Pr(0,1,3)*dc[3]);
  dcP[ 3] = (Pi(0,1,0)*dc[0]+Pi(0,1,1)*dc[1]+Pi(0,1,2)*dc[2]+Pi(0,1,3)*dc[3]);
  dcP[ 4] = (Pr(0,2,0)*dc[0]+Pr(0,2,1)*dc[1]+Pr(0,2,2)*dc[2]+Pr(0,2,3)*dc[3]);
  dcP[ 5] = (Pi(0,2,0)*dc[0]+Pi(0,2,1)*dc[1]+Pi(0,2,2)*dc[2]+Pi(0,2,3)*dc[3]);
  dcP[ 6] = (Pr(0,3,0)*dc[0]+Pr(0,3,1)*dc[1]+Pr(0,3,2)*dc[2]+Pr(0,3,3)*dc[3]);
  dcP[ 7] = (Pi(0,3,0)*dc[0]+Pi(0,3,1)*dc[1]+Pi(0,3,2)*dc[2]+Pi(0,3,3)*dc[3]);
  dcP[ 8] = (Pr(1,0,0)*dc[0]+Pr(1,0,1)*dc[1]+Pr(1,0,2)*dc[2]+Pr(1,0,3)*dc[3]);
  dcP[ 9] = (Pi(1,0,0)*dc[0]+Pi(1,0,1)*dc[1]+Pi(1,0,2)*dc[2]+Pi(1,0,3)*dc[3]);
  dcP[10] = (Pr(1,1,0)*dc[0]+Pr(1,1,1)*dc[1]+Pr(1,1,2)*dc[2]+Pr(1,1,3)*dc[3]);
  dcP[11] = (Pi(1,1,0)*dc[0]+Pi(1,1,1)*dc[1]+Pi(1,1,2)*dc[2]+Pi(1,1,3)*dc[3]);
  dcP[12] = (Pr(1,2,0)*dc[0]+Pr(1,2,1)*dc[1]+Pr(1,2,2)*dc[2]+Pr(1,2,3)*dc[3]);
  dcP[13] = (Pi(1,2,0)*dc[0]+Pi(1,2,1)*dc[1]+Pi(1,2,2)*dc[2]+Pi(1,2,3)*dc[3]);
  dcP[14] = (Pr(1,3,0)*dc[0]+Pr(1,3,1)*dc[1]+Pr(1,3,2)*dc[2]+Pr(1,3,3)*dc[3]);
  dcP[15] = (Pi(1,3,0)*dc[0]+Pi(1,3,1)*dc[1]+Pi(1,3,2)*dc[2]+Pi(1,3,3)*dc[3]);
  dcP[16] = (Pr(2,0,0)*dc[0]+Pr(2,0,1)*dc[1]+Pr(2,0,2)*dc[2]+Pr(2,0,3)*dc[3]);
  dcP[17] = (Pi(2,0,0)*dc[0]+Pi(2,0,1)*dc[1]+Pi(2,0,2)*dc[2]+Pi(2,0,3)*dc[3]);
  dcP[18] = (Pr(2,1,0)*dc[0]+Pr(2,1,1)*dc[1]+Pr(2,1,2)*dc[2]+Pr(2,1,3)*dc[3]);
  dcP[19] = (Pi(2,1,0)*dc[0]+Pi(2,1,1)*dc[1]+Pi(2,1,2)*dc[2]+Pi(2,1,3)*dc[3]);
  dcP[20] = (Pr(2,2,0)*dc[0]+Pr(2,2,1)*dc[1]+Pr(2,2,2)*dc[2]+Pr(2,2,3)*dc[3]);
  dcP[21] = (Pi(2,2,0)*dc[0]+Pi(2,2,1)*dc[1]+Pi(2,2,2)*dc[2]+Pi(2,2,3)*dc[3]);
  dcP[22] = (Pr(2,3,0)*dc[0]+Pr(2,3,1)*dc[1]+Pr(2,3,2)*dc[2]+Pr(2,3,3)*dc[3]);
  dcP[23] = (Pi(2,3,0)*dc[0]+Pi(2,3,1)*dc[1]+Pi(2,3,2)*dc[2]+Pi(2,3,3)*dc[3]);
  dcP[24] = (Pr(3,0,0)*dc[0]+Pr(3,0,1)*dc[1]+Pr(3,0,2)*dc[2]+Pr(3,0,3)*dc[3]);
  dcP[25] = (Pi(3,0,0)*dc[0]+Pi(3,0,1)*dc[1]+Pi(3,0,2)*dc[2]+Pi(3,0,3)*dc[3]);
  dcP[26] = (Pr(3,1,0)*dc[0]+Pr(3,1,1)*dc[1]+Pr(3,1,2)*dc[2]+Pr(3,1,3)*dc[3]);
  dcP[27] = (Pi(3,1,0)*dc[0]+Pi(3,1,1)*dc[1]+Pi(3,1,2)*dc[2]+Pi(3,1,3)*dc[3]);
  dcP[28] = (Pr(3,2,0)*dc[0]+Pr(3,2,1)*dc[1]+Pr(3,2,2)*dc[2]+Pr(3,2,3)*dc[3]);
  dcP[29] = (Pi(3,2,0)*dc[0]+Pi(3,2,1)*dc[1]+Pi(3,2,2)*dc[2]+Pi(3,2,3)*dc[3]);
  dcP[30] = (Pr(3,3,0)*dc[0]+Pr(3,3,1)*dc[1]+Pr(3,3,2)*dc[2]+Pr(3,3,3)*dc[3]);
  dcP[31] = (Pi(3,3,0)*dc[0]+Pi(3,3,1)*dc[1]+Pi(3,3,2)*dc[2]+Pi(3,3,3)*dc[3]);
  bcP[0] = ( b[0]*cP[ 0] +  b[1]*cP[ 2] +  b[2]*cP[ 4] +  b[3]*cP[ 6]);
  bcP[1] = ( b[0]*cP[ 1] +  b[1]*cP[ 3] +  b[2]*cP[ 5] +  b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] +  b[1]*cP[10] +  b[2]*cP[12] +  b[3]*cP[14]);
  bcP[3] = ( b[0]*cP[ 9] +  b[1]*cP[11] +  b[2]*cP[13] +  b[3]*cP[15]);
  bcP[4] = ( b[0]*cP[16] +  b[1]*cP[18] +  b[2]*cP[20] +  b[3]*cP[22]);
  bcP[5] = ( b[0]*cP[17] +  b[1]*cP[19] +  b[2]*cP[21] +  b[3]*cP[23]);
  bcP[6] = ( b[0]*cP[24] +  b[1]*cP[26] +  b[2]*cP[28] +  b[3]*cP[30]);
  bcP[7] = ( b[0]*cP[25] +  b[1]*cP[27] +  b[2]*cP[29] +  b[3]*cP[31]);
  dbcP[0] = (db[0]*cP[ 0] + db[1]*cP[ 2] + db[2]*cP[ 4] + db[3]*cP[ 6]);
  dbcP[1] = (db[0]*cP[ 1] + db[1]*cP[ 3] + db[2]*cP[ 5] + db[3]*cP[ 7]);
  dbcP[2] = (db[0]*cP[ 8] + db[1]*cP[10] + db[2]*cP[12] + db[3]*cP[14]);
  dbcP[3] = (db[0]*cP[ 9] + db[1]*cP[11] + db[2]*cP[13] + db[3]*cP[15]);
  dbcP[4] = (db[0]*cP[16] + db[1]*cP[18] + db[2]*cP[20] + db[3]*cP[22]);
  dbcP[5] = (db[0]*cP[17] + db[1]*cP[19] + db[2]*cP[21] + db[3]*cP[23]);
  dbcP[6] = (db[0]*cP[24] + db[1]*cP[26] + db[2]*cP[28] + db[3]*cP[30]);
  dbcP[7] = (db[0]*cP[25] + db[1]*cP[27] + db[2]*cP[29] + db[3]*cP[31]);
  bdcP[0] = (b[0]*dcP[ 0] + b[1]*dcP[ 2] + b[2]*dcP[ 4] + b[3]*dcP[ 6]);
  bdcP[1] = (b[0]*dcP[ 1] + b[1]*dcP[ 3] + b[2]*dcP[ 5] + b[3]*dcP[ 7]);
  bdcP[2] = (b[0]*dcP[ 8] + b[1]*dcP[10] + b[2]*dcP[12] + b[3]*dcP[14]);
  bdcP[3] = (b[0]*dcP[ 9] + b[1]*dcP[11] + b[2]*dcP[13] + b[3]*dcP[15]);
  bdcP[4] = (b[0]*dcP[16] + b[1]*dcP[18] + b[2]*dcP[20] + b[3]*dcP[22]);
  bdcP[5] = (b[0]*dcP[17] + b[1]*dcP[19] + b[2]*dcP[21] + b[3]*dcP[23]);
  bdcP[6] = (b[0]*dcP[24] + b[1]*dcP[26] + b[2]*dcP[28] + b[3]*dcP[30]);
  bdcP[7] = (b[0]*dcP[25] + b[1]*dcP[27] + b[2]*dcP[29] + b[3]*dcP[31]);
  val[0]    = ( a[0]*bcP[0] +  a[1]*bcP[2] +  a[2]*bcP[4] +  a[3]*bcP[6]);
  val[1]    = ( a[0]*bcP[1] +  a[1]*bcP[3] +  a[2]*bcP[5] +  a[3]*bcP[7]);
  grad[0] = spline->x_grid.delta_inv *
            (da[0]*bcP[0] + da[1]*bcP[2] + da[2]*bcP[4] + da[3]*bcP[6]);
  grad[1] = spline->x_grid.delta_inv *
            (da[0]*bcP[1] + da[1]*bcP[3] + da[2]*bcP[5] + da[3]*bcP[7]);
  grad[2] = spline->y_grid.delta_inv *
            (a[0]*dbcP[0] + a[1]*dbcP[2] + a[2]*dbcP[4] + a[3]*dbcP[6]);
  grad[3] = spline->y_grid.delta_inv *
            (a[0]*dbcP[1] + a[1]*dbcP[3] + a[2]*dbcP[5] + a[3]*dbcP[7]);
  grad[4] = spline->z_grid.delta_inv *
            (a[0]*bdcP[0] + a[1]*bdcP[2] + a[2]*bdcP[4] + a[3]*bdcP[6]);
  grad[5] = spline->z_grid.delta_inv *
            (a[5]*bdcP[1] + a[1]*bdcP[3] + a[2]*bdcP[5] + a[3]*bdcP[7]);
#undef Pr
#undef Pi
}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_c_vgl (UBspline_3d_c * restrict spline,
                        double x, double y, double z,
                        float* restrict val, float* restrict grad, float* restrict lapl)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modff (uy, &iparty);
  int iy = (int) iparty;
  tz = modff (uz, &ipartz);
  int iz = (int) ipartz;
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4],
        d2a[4], d2b[4], d2c[4], cP[16], dcP[16], bcP[4], dbcP[4], d2bcP[4], bdcP[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;
  tpz[1] = tz*tz;
  tpz[2] = tz;
  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);
  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);
  int xs = spline->x_ctride;
  int ys = spline->y_ctride;
#define Pr(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))]
#define Pi(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))+1]
  cP[ 0] = (P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3]);
  cP[ 1] = (P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3]);
  cP[ 2] = (P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3]);
  cP[ 3] = (P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]);
  cP[ 4] = (P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3]);
  cP[ 5] = (P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3]);
  cP[ 6] = (P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3]);
  cP[ 7] = (P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]);
  cP[ 8] = (P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3]);
  cP[ 9] = (P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3]);
  cP[10] = (P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3]);
  cP[11] = (P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]);
  cP[12] = (P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3]);
  cP[13] = (P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3]);
  cP[14] = (P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3]);
  cP[15] = (P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3]);
  dcP[ 0] = (P(0,0,0)*dc[0]+P(0,0,1)*dc[1]+P(0,0,2)*dc[2]+P(0,0,3)*dc[3]);
  dcP[ 1] = (P(0,1,0)*dc[0]+P(0,1,1)*dc[1]+P(0,1,2)*dc[2]+P(0,1,3)*dc[3]);
  dcP[ 2] = (P(0,2,0)*dc[0]+P(0,2,1)*dc[1]+P(0,2,2)*dc[2]+P(0,2,3)*dc[3]);
  dcP[ 3] = (P(0,3,0)*dc[0]+P(0,3,1)*dc[1]+P(0,3,2)*dc[2]+P(0,3,3)*dc[3]);
  dcP[ 4] = (P(1,0,0)*dc[0]+P(1,0,1)*dc[1]+P(1,0,2)*dc[2]+P(1,0,3)*dc[3]);
  dcP[ 5] = (P(1,1,0)*dc[0]+P(1,1,1)*dc[1]+P(1,1,2)*dc[2]+P(1,1,3)*dc[3]);
  dcP[ 6] = (P(1,2,0)*dc[0]+P(1,2,1)*dc[1]+P(1,2,2)*dc[2]+P(1,2,3)*dc[3]);
  dcP[ 7] = (P(1,3,0)*dc[0]+P(1,3,1)*dc[1]+P(1,3,2)*dc[2]+P(1,3,3)*dc[3]);
  dcP[ 8] = (P(2,0,0)*dc[0]+P(2,0,1)*dc[1]+P(2,0,2)*dc[2]+P(2,0,3)*dc[3]);
  dcP[ 9] = (P(2,1,0)*dc[0]+P(2,1,1)*dc[1]+P(2,1,2)*dc[2]+P(2,1,3)*dc[3]);
  dcP[10] = (P(2,2,0)*dc[0]+P(2,2,1)*dc[1]+P(2,2,2)*dc[2]+P(2,2,3)*dc[3]);
  dcP[11] = (P(2,3,0)*dc[0]+P(2,3,1)*dc[1]+P(2,3,2)*dc[2]+P(2,3,3)*dc[3]);
  dcP[12] = (P(3,0,0)*dc[0]+P(3,0,1)*dc[1]+P(3,0,2)*dc[2]+P(3,0,3)*dc[3]);
  dcP[13] = (P(3,1,0)*dc[0]+P(3,1,1)*dc[1]+P(3,1,2)*dc[2]+P(3,1,3)*dc[3]);
  dcP[14] = (P(3,2,0)*dc[0]+P(3,2,1)*dc[1]+P(3,2,2)*dc[2]+P(3,2,3)*dc[3]);
  dcP[15] = (P(3,3,0)*dc[0]+P(3,3,1)*dc[1]+P(3,3,2)*dc[2]+P(3,3,3)*dc[3]);
  bcP[0] = ( b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3]);
  bcP[1] = ( b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11]);
  bcP[3] = ( b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15]);
  dbcP[0] = ( db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]);
  dbcP[1] = ( db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]);
  dbcP[2] = ( db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]);
  dbcP[3] = ( db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]);
  bdcP[0] = ( b[0]*dcP[ 0] + b[1]*dcP[ 1] + b[2]*dcP[ 2] + b[3]*dcP[ 3]);
  bdcP[1] = ( b[0]*dcP[ 4] + b[1]*dcP[ 5] + b[2]*dcP[ 6] + b[3]*dcP[ 7]);
  bdcP[2] = ( b[0]*dcP[ 8] + b[1]*dcP[ 9] + b[2]*dcP[10] + b[3]*dcP[11]);
  bdcP[3] = ( b[0]*dcP[12] + b[1]*dcP[13] + b[2]*dcP[14] + b[3]*dcP[15]);
  d2bcP[0] = ( d2b[0]*cP[ 0] + d2b[1]*cP[ 1] + d2b[2]*cP[ 2] + d2b[3]*cP[ 3]);
  d2bcP[1] = ( d2b[0]*cP[ 4] + d2b[1]*cP[ 5] + d2b[2]*cP[ 6] + d2b[3]*cP[ 7]);
  d2bcP[2] = ( d2b[0]*cP[ 8] + d2b[1]*cP[ 9] + d2b[2]*cP[10] + d2b[3]*cP[11]);
  d2bcP[3] = ( d2b[0]*cP[12] + d2b[1]*cP[13] + d2b[2]*cP[14] + d2b[3]*cP[15]);
  *val    =
    ( a[0]*bcP[0] +  a[1]*bcP[1] +  a[2]*bcP[2] +  a[3]*bcP[3]);
  grad[0] = spline->x_grid.delta_inv *
            (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = spline->y_grid.delta_inv *
            (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = spline->z_grid.delta_inv *
            (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);
  *lapl =
    spline->x_grid.delta_inv * spline->x_grid.delta_inv *
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3])
    + spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]) +
    + spline->z_grid.delta_inv * spline->z_grid.delta_inv *
    (a[0]*(b[0]*(P(0,0,0)*d2c[0]+P(0,0,1)*d2c[1]+P(0,0,2)*d2c[2]+P(0,0,3)*d2c[3])+
           b[1]*(P(0,1,0)*d2c[0]+P(0,1,1)*d2c[1]+P(0,1,2)*d2c[2]+P(0,1,3)*d2c[3])+
           b[2]*(P(0,2,0)*d2c[0]+P(0,2,1)*d2c[1]+P(0,2,2)*d2c[2]+P(0,2,3)*d2c[3])+
           b[3]*(P(0,3,0)*d2c[0]+P(0,3,1)*d2c[1]+P(0,3,2)*d2c[2]+P(0,3,3)*d2c[3]))+
     a[1]*(b[0]*(P(1,0,0)*d2c[0]+P(1,0,1)*d2c[1]+P(1,0,2)*d2c[2]+P(1,0,3)*d2c[3])+
           b[1]*(P(1,1,0)*d2c[0]+P(1,1,1)*d2c[1]+P(1,1,2)*d2c[2]+P(1,1,3)*d2c[3])+
           b[2]*(P(1,2,0)*d2c[0]+P(1,2,1)*d2c[1]+P(1,2,2)*d2c[2]+P(1,2,3)*d2c[3])+
           b[3]*(P(1,3,0)*d2c[0]+P(1,3,1)*d2c[1]+P(1,3,2)*d2c[2]+P(1,3,3)*d2c[3]))+
     a[2]*(b[0]*(P(2,0,0)*d2c[0]+P(2,0,1)*d2c[1]+P(2,0,2)*d2c[2]+P(2,0,3)*d2c[3])+
           b[1]*(P(2,1,0)*d2c[0]+P(2,1,1)*d2c[1]+P(2,1,2)*d2c[2]+P(2,1,3)*d2c[3])+
           b[2]*(P(2,2,0)*d2c[0]+P(2,2,1)*d2c[1]+P(2,2,2)*d2c[2]+P(2,2,3)*d2c[3])+
           b[3]*(P(2,3,0)*d2c[0]+P(2,3,1)*d2c[1]+P(2,3,2)*d2c[2]+P(2,3,3)*d2c[3]))+
     a[3]*(b[0]*(P(3,0,0)*d2c[0]+P(3,0,1)*d2c[1]+P(3,0,2)*d2c[2]+P(3,0,3)*d2c[3])+
           b[1]*(P(3,1,0)*d2c[0]+P(3,1,1)*d2c[1]+P(3,1,2)*d2c[2]+P(3,1,3)*d2c[3])+
           b[2]*(P(3,2,0)*d2c[0]+P(3,2,1)*d2c[1]+P(3,2,2)*d2c[2]+P(3,2,3)*d2c[3])+
           b[3]*(P(3,3,0)*d2c[0]+P(3,3,1)*d2c[1]+P(3,3,2)*d2c[2]+P(3,3,3)*d2c[3])));
#undef Pr
#undef Pi
}





/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_c_vgh (UBspline_3d_c * restrict spline,
                        double x, double y, double z,
                        float* restrict val, float* restrict grad, float* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modff (uy, &iparty);
  int iy = (int) iparty;
  tz = modff (uz, &ipartz);
  int iz = (int) ipartz;
//   if ((ix >= spline->x_grid.num))    x = spline->x_grid.num;
//   if ((ix < 0))                      x = 0;
//   if ((iy >= spline->y_grid.num))    y = spline->y_grid.num;
//   if ((iy < 0))                      y = 0;
//   if ((iz >= spline->z_grid.num))    z = spline->z_grid.num;
//   if ((iz < 0))                      z = 0;
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4],
        d2a[4], d2b[4], d2c[4], cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
        d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;
  tpz[1] = tz*tz;
  tpz[2] = tz;
  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;
  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);
  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);
  int xs = spline->x_ctride;
  int ys = spline->y_ctride;
#define Pr(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))]
#define Pi(i,j,k) coefs[2*((ix+(i))*xs+(iy+(j))*ys+(iz+(k)))+1]
  cP[ 0] = (P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3]);
  cP[ 1] = (P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3]);
  cP[ 2] = (P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3]);
  cP[ 3] = (P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]);
  cP[ 4] = (P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3]);
  cP[ 5] = (P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3]);
  cP[ 6] = (P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3]);
  cP[ 7] = (P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]);
  cP[ 8] = (P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3]);
  cP[ 9] = (P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3]);
  cP[10] = (P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3]);
  cP[11] = (P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]);
  cP[12] = (P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3]);
  cP[13] = (P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3]);
  cP[14] = (P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3]);
  cP[15] = (P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3]);
  dcP[ 0] = (P(0,0,0)*dc[0]+P(0,0,1)*dc[1]+P(0,0,2)*dc[2]+P(0,0,3)*dc[3]);
  dcP[ 1] = (P(0,1,0)*dc[0]+P(0,1,1)*dc[1]+P(0,1,2)*dc[2]+P(0,1,3)*dc[3]);
  dcP[ 2] = (P(0,2,0)*dc[0]+P(0,2,1)*dc[1]+P(0,2,2)*dc[2]+P(0,2,3)*dc[3]);
  dcP[ 3] = (P(0,3,0)*dc[0]+P(0,3,1)*dc[1]+P(0,3,2)*dc[2]+P(0,3,3)*dc[3]);
  dcP[ 4] = (P(1,0,0)*dc[0]+P(1,0,1)*dc[1]+P(1,0,2)*dc[2]+P(1,0,3)*dc[3]);
  dcP[ 5] = (P(1,1,0)*dc[0]+P(1,1,1)*dc[1]+P(1,1,2)*dc[2]+P(1,1,3)*dc[3]);
  dcP[ 6] = (P(1,2,0)*dc[0]+P(1,2,1)*dc[1]+P(1,2,2)*dc[2]+P(1,2,3)*dc[3]);
  dcP[ 7] = (P(1,3,0)*dc[0]+P(1,3,1)*dc[1]+P(1,3,2)*dc[2]+P(1,3,3)*dc[3]);
  dcP[ 8] = (P(2,0,0)*dc[0]+P(2,0,1)*dc[1]+P(2,0,2)*dc[2]+P(2,0,3)*dc[3]);
  dcP[ 9] = (P(2,1,0)*dc[0]+P(2,1,1)*dc[1]+P(2,1,2)*dc[2]+P(2,1,3)*dc[3]);
  dcP[10] = (P(2,2,0)*dc[0]+P(2,2,1)*dc[1]+P(2,2,2)*dc[2]+P(2,2,3)*dc[3]);
  dcP[11] = (P(2,3,0)*dc[0]+P(2,3,1)*dc[1]+P(2,3,2)*dc[2]+P(2,3,3)*dc[3]);
  dcP[12] = (P(3,0,0)*dc[0]+P(3,0,1)*dc[1]+P(3,0,2)*dc[2]+P(3,0,3)*dc[3]);
  dcP[13] = (P(3,1,0)*dc[0]+P(3,1,1)*dc[1]+P(3,1,2)*dc[2]+P(3,1,3)*dc[3]);
  dcP[14] = (P(3,2,0)*dc[0]+P(3,2,1)*dc[1]+P(3,2,2)*dc[2]+P(3,2,3)*dc[3]);
  dcP[15] = (P(3,3,0)*dc[0]+P(3,3,1)*dc[1]+P(3,3,2)*dc[2]+P(3,3,3)*dc[3]);
  d2cP[ 0] = (P(0,0,0)*d2c[0]+P(0,0,1)*d2c[1]+P(0,0,2)*d2c[2]+P(0,0,3)*d2c[3]);
  d2cP[ 1] = (P(0,1,0)*d2c[0]+P(0,1,1)*d2c[1]+P(0,1,2)*d2c[2]+P(0,1,3)*d2c[3]);
  d2cP[ 2] = (P(0,2,0)*d2c[0]+P(0,2,1)*d2c[1]+P(0,2,2)*d2c[2]+P(0,2,3)*d2c[3]);
  d2cP[ 3] = (P(0,3,0)*d2c[0]+P(0,3,1)*d2c[1]+P(0,3,2)*d2c[2]+P(0,3,3)*d2c[3]);
  d2cP[ 4] = (P(1,0,0)*d2c[0]+P(1,0,1)*d2c[1]+P(1,0,2)*d2c[2]+P(1,0,3)*d2c[3]);
  d2cP[ 5] = (P(1,1,0)*d2c[0]+P(1,1,1)*d2c[1]+P(1,1,2)*d2c[2]+P(1,1,3)*d2c[3]);
  d2cP[ 6] = (P(1,2,0)*d2c[0]+P(1,2,1)*d2c[1]+P(1,2,2)*d2c[2]+P(1,2,3)*d2c[3]);
  d2cP[ 7] = (P(1,3,0)*d2c[0]+P(1,3,1)*d2c[1]+P(1,3,2)*d2c[2]+P(1,3,3)*d2c[3]);
  d2cP[ 8] = (P(2,0,0)*d2c[0]+P(2,0,1)*d2c[1]+P(2,0,2)*d2c[2]+P(2,0,3)*d2c[3]);
  d2cP[ 9] = (P(2,1,0)*d2c[0]+P(2,1,1)*d2c[1]+P(2,1,2)*d2c[2]+P(2,1,3)*d2c[3]);
  d2cP[10] = (P(2,2,0)*d2c[0]+P(2,2,1)*d2c[1]+P(2,2,2)*d2c[2]+P(2,2,3)*d2c[3]);
  d2cP[11] = (P(2,3,0)*d2c[0]+P(2,3,1)*d2c[1]+P(2,3,2)*d2c[2]+P(2,3,3)*d2c[3]);
  d2cP[12] = (P(3,0,0)*d2c[0]+P(3,0,1)*d2c[1]+P(3,0,2)*d2c[2]+P(3,0,3)*d2c[3]);
  d2cP[13] = (P(3,1,0)*d2c[0]+P(3,1,1)*d2c[1]+P(3,1,2)*d2c[2]+P(3,1,3)*d2c[3]);
  d2cP[14] = (P(3,2,0)*d2c[0]+P(3,2,1)*d2c[1]+P(3,2,2)*d2c[2]+P(3,2,3)*d2c[3]);
  d2cP[15] = (P(3,3,0)*d2c[0]+P(3,3,1)*d2c[1]+P(3,3,2)*d2c[2]+P(3,3,3)*d2c[3]);
  bcP[0] = ( b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3]);
  bcP[1] = ( b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11]);
  bcP[3] = ( b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15]);
  dbcP[0] = ( db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]);
  dbcP[1] = ( db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]);
  dbcP[2] = ( db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]);
  dbcP[3] = ( db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]);
  bdcP[0] = ( b[0]*dcP[ 0] + b[1]*dcP[ 1] + b[2]*dcP[ 2] + b[3]*dcP[ 3]);
  bdcP[1] = ( b[0]*dcP[ 4] + b[1]*dcP[ 5] + b[2]*dcP[ 6] + b[3]*dcP[ 7]);
  bdcP[2] = ( b[0]*dcP[ 8] + b[1]*dcP[ 9] + b[2]*dcP[10] + b[3]*dcP[11]);
  bdcP[3] = ( b[0]*dcP[12] + b[1]*dcP[13] + b[2]*dcP[14] + b[3]*dcP[15]);
  bd2cP[0] = ( b[0]*d2cP[ 0] + b[1]*d2cP[ 1] + b[2]*d2cP[ 2] + b[3]*d2cP[ 3]);
  bd2cP[1] = ( b[0]*d2cP[ 4] + b[1]*d2cP[ 5] + b[2]*d2cP[ 6] + b[3]*d2cP[ 7]);
  bd2cP[2] = ( b[0]*d2cP[ 8] + b[1]*d2cP[ 9] + b[2]*d2cP[10] + b[3]*d2cP[11]);
  bd2cP[3] = ( b[0]*d2cP[12] + b[1]*d2cP[13] + b[2]*d2cP[14] + b[3]*d2cP[15]);
  d2bcP[0] = ( d2b[0]*cP[ 0] + d2b[1]*cP[ 1] + d2b[2]*cP[ 2] + d2b[3]*cP[ 3]);
  d2bcP[1] = ( d2b[0]*cP[ 4] + d2b[1]*cP[ 5] + d2b[2]*cP[ 6] + d2b[3]*cP[ 7]);
  d2bcP[2] = ( d2b[0]*cP[ 8] + d2b[1]*cP[ 9] + d2b[2]*cP[10] + d2b[3]*cP[11]);
  d2bcP[3] = ( d2b[0]*cP[12] + d2b[1]*cP[13] + d2b[2]*cP[14] + d2b[3]*cP[15]);
  dbdcP[0] = ( db[0]*dcP[ 0] + db[1]*dcP[ 1] + db[2]*dcP[ 2] + db[3]*dcP[ 3]);
  dbdcP[1] = ( db[0]*dcP[ 4] + db[1]*dcP[ 5] + db[2]*dcP[ 6] + db[3]*dcP[ 7]);
  dbdcP[2] = ( db[0]*dcP[ 8] + db[1]*dcP[ 9] + db[2]*dcP[10] + db[3]*dcP[11]);
  dbdcP[3] = ( db[0]*dcP[12] + db[1]*dcP[13] + db[2]*dcP[14] + db[3]*dcP[15]);
  *val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  grad[0] = spline->x_grid.delta_inv *
            (da[0] *bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = spline->y_grid.delta_inv *
            (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = spline->z_grid.delta_inv *
            (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);
  // d2x
  hess[0] = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
            (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3]);
  // dx dy
  hess[1] = spline->x_grid.delta_inv * spline->y_grid.delta_inv *
            (da[0]*dbcP[0] + da[1]*dbcP[1] + da[1]*dbcP[1] + da[1]*dbcP[1]);
  hess[3] = hess[1];
  // dx dz;
  hess[2] = spline->x_grid.delta_inv * spline->z_grid.delta_inv *
            (da[0]*bdcP[0] + da[1]*bdcP[1] + da[1]*bdcP[1] + da[1]*bdcP[1]);
  hess[6] = hess[2];
  // d2y
  hess[4] = spline->y_grid.delta_inv * spline->y_grid.delta_inv *
            (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]);
  // dy dz
  hess[5] = spline->y_grid.delta_inv * spline->z_grid.delta_inv *
            (a[0]*dbdcP[0] + a[1]*dbdcP[1] + a[2]*dbdcP[2] + a[3]*dbdcP[3]);
  hess[7] = hess[5];
  // d2z
  hess[8] = spline->z_grid.delta_inv * spline->z_grid.delta_inv *
            (a[0]*bd2cP[0] + a[1]*bd2cP[1] + a[2]*bd2cP[2] + a[3]*bd2cP[3]);
#undef Pr
#undef Pi
}

#endif
