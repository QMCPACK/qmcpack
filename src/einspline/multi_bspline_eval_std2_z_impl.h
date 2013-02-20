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

#ifndef MULTI_BSPLINE_EVAL_STD2_Z_IMPL_H
#define MULTI_BSPLINE_EVAL_STD2_Z_IMPL_H

#include <math.h>
#include <stdio.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"

extern const double* restrict   Ad;
extern const double* restrict  dAd;
extern const double* restrict d2Ad;
extern const double* restrict d3Ad;

/************************************************************/
/* 1D double-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_z (const multi_UBspline_1d_z *spline,
			  double x,
			  complex_double* restrict vals)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;
  
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  for (int n=0; n<spline->num_splines; n++) 
    vals[n]  = 0.0;

  for (int i=0; i<4; i++) {
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++) 
      vals[n]  +=   a[i] * coefs[n];
  }
}



void
eval_multi_UBspline_1d_z_vg (const multi_UBspline_1d_z *spline,
			     double x,
			     complex_double* restrict vals,
			     complex_double* restrict grads)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

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
  
  double dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) 
    grads[n] *= dxInv;
}


void
eval_multi_UBspline_1d_z_vgl (const multi_UBspline_1d_z *spline,
			      double x,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict lapl)	  
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

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

  double dxInv = spline->x_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[n] *= dxInv;
    lapl [n] *= dxInv*dxInv;
  }
}


void
eval_multi_UBspline_1d_z_vgh (const multi_UBspline_1d_z *spline,
			      double x,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict hess)
{
  eval_multi_UBspline_1d_z_vgl (spline, x, vals, grads, hess);
}


/************************************************************/
/* 2D double-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_z (const multi_UBspline_2d_z *spline,
			  double x, double y,
			  complex_double* restrict vals)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  double tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      double prefactor = a[i]*b[j];
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) 
	vals[n] += prefactor*coefs[n];
    }
}


void
eval_multi_UBspline_2d_z_vg (const multi_UBspline_2d_z *spline,
			     double x, double y,
			     complex_double* restrict vals,
			     complex_double* restrict grads)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = grads[2*n+2] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      double ab = a[i]*b[j];
      double dab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals [n]     +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
      }
    }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
  }
}

void
eval_multi_UBspline_2d_z_vgl (const multi_UBspline_2d_z *spline,
			      double x, double y,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict lapl)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

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

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  //complex_double lapl2[2*spline->num_splines];
  complex_double* restrict lapl2 = spline->lapl2;
  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    lapl2[2*n+0] = lapl2[2*n+1] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      double ab = a[i]*b[j];
      double dab[2], d2ab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i]*  b[j];
      d2ab[1] =   a[i]*d2b[j];
      
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals[n]      +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
	lapl2[2*n+0] += d2ab[0]*coefs[n];
	lapl2[2*n+1] += d2ab[1]*coefs[n];
      }
    }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    lapl2[2*n+0] *= dxInv*dxInv;
    lapl2[2*n+1] *= dyInv*dyInv;
    lapl[n] = lapl2[2*n+0] + lapl2[2*n+1];
  }
}




void
eval_multi_UBspline_2d_z_vgh (const multi_UBspline_2d_z *spline,
			      double x, double y,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict hess)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

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

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    for (int i=0; i<4; i++)
      hess[4*n+i] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      double ab = a[i]*b[j];
      double dab[2], d2ab[3];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i] *   b[j];
      d2ab[1] =  da[i] *  db[j];
      d2ab[2] =   a[i] * d2b[j];

      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals[n]      +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
	hess [4*n+0] += d2ab[0]*coefs[n];
	hess [4*n+1] += d2ab[1]*coefs[n];
	hess [4*n+3] += d2ab[2]*coefs[n];
      }
    }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    hess[4*n+0] *= dxInv*dxInv;
    hess[4*n+1] *= dxInv*dyInv;
    hess[4*n+3] *= dyInv*dyInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess[4*n+2] = hess[4*n+1];
  }
}



/************************************************************/
/* 3D double-precision, complex evaulation functions        */
/************************************************************/


void 
eval_multi_UBspline_3d_z(const multi_UBspline_3d_z *spline, 
                         double x, double y, double z, 
                         complex_double* restrict vals)
{

    double ux, uy, uz, ipartx, iparty, ipartz, tx, ty, tz, a[4], b[4],c[4], d[64], s;
    complex_double *mod_coefs[64];
    intptr_t xs, ys, zs;
    int i, j, k, n, M;
    double *v, *p;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;
    
    tx = modf (ux, &ipartx);  int ix = (int) ipartx;
    ty = modf (uy, &iparty);  int iy = (int) iparty;
    tz = modf (uz, &ipartz);  int iz = (int) ipartz;

    a[0] = ( ( Ad[0]  * tx + Ad[1] ) * tx + Ad[2] ) * tx + Ad[3]; 
    a[1] = ( ( Ad[4]  * tx + Ad[5] ) * tx + Ad[6] ) * tx + Ad[7]; 
    a[2] = ( ( Ad[8]  * tx + Ad[9] ) * tx + Ad[10] ) * tx + Ad[11]; 
    a[3] = ( ( Ad[12] * tx + Ad[13] ) * tx + Ad[14] ) * tx + Ad[15]; 
    
    b[0] = ( ( Ad[0]  * ty + Ad[1] ) * ty + Ad[2] ) * ty + Ad[3]; 
    b[1] = ( ( Ad[4]  * ty + Ad[5] ) * ty + Ad[6] ) * ty + Ad[7]; 
    b[2] = ( ( Ad[8]  * ty + Ad[9] ) * ty + Ad[10] ) * ty + Ad[11]; 
    b[3] = ( ( Ad[12] * ty + Ad[13] ) * ty + Ad[14] ) * ty + Ad[15]; 
    
    c[0] = ( ( Ad[0]  * tz + Ad[1] ) * tz + Ad[2] ) * tz + Ad[3]; 
    c[1] = ( ( Ad[4]  * tz + Ad[5] ) * tz + Ad[6] ) * tz + Ad[7]; 
    c[2] = ( ( Ad[8]  * tz + Ad[9] ) * tz + Ad[10] ) * tz + Ad[11]; 
    c[3] = ( ( Ad[12] * tz + Ad[13] ) * tz + Ad[14] ) * tz + Ad[15]; 

    xs = spline->x_stride;
    ys = spline->y_stride;
    zs = spline->z_stride;
    n = 0;
    for ( i=0; i<4; i++)
      for ( j=0; j<4; j++) 
        for ( k=0; k<4; k++) {
	  d[n] = a[i]*b[j]*c[k];
          mod_coefs[n] = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	  n = n + 1;
        }
      

    M = 2 * ( spline -> num_splines );
    v = (double *)&vals[0];
    int rem = M % 8;    
    for ( n = 0; n < M; n++ ) v[n] = 0.0;

    for ( i = 0; i < 64; i++ )
    {  
       s = d[i];
       p = (double *)mod_coefs[i];
       
       for ( n = 0; n < M - rem ; n = n + 8 ) 
       {
          double a0, a1, a2, a3, a4, a5, a6, a7;
          double b0, b1, b2, b3, b4, b5, b6, b7;

          a0 = v[n];
          a1 = v[n+1];
          a2 = v[n+2];
          a3 = v[n+3];
          a4 = v[n+4];
          a5 = v[n+5];
          a6 = v[n+6];
          a7 = v[n+7];
          
          b0 = p[n];
          b1 = p[n+1];
          b2 = p[n+2];
          b3 = p[n+3];
          b4 = p[n+4];
          b5 = p[n+5];
          b6 = p[n+6];
          b7 = p[n+7];
          
          a0 = a0 + s * b0;
          a1 = a1 + s * b1;
          a2 = a2 + s * b2;
          a3 = a3 + s * b3;
          a4 = a4 + s * b4;
          a5 = a5 + s * b5;
          a6 = a6 + s * b6;
          a7 = a7 + s * b7;
 
          v[n  ] = a0;
          v[n+1] = a1;
          v[n+2] = a2;
          v[n+3] = a3;
          v[n+4] = a4;
          v[n+5] = a5;
          v[n+6] = a6;
          v[n+7] = a7;
      }

       for ( k = n; k < M; k++ ){
           v[k] = v[k] + s * p[k];
       }

    } 
}

void
eval_multi_UBspline_3d_z_vgh (const multi_UBspline_3d_z *spline,
			      double x, double y, double z,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict hess)	  
{

    double ux, uy, uz, prefactor, ipartx, iparty, ipartz, tx, ty, tz;
    double dxInv,dxInvdxInv,dyInv,dyInvdyInv,dzInv,dzInvdzInv,dxInvdyInv,dxInvdzInv,dyInvdzInv;
    double  a[4], b[4],c[4], da[4], db[4],dc[4], d2a[4], d2b[4],d2c[4]; 
    double dVal[64],dGrad0[64],dGrad1[64],dGrad2[64],dHess0[64],dHess1[64],dHess2[64],dHess4[64],dHess5[64],dHess8[64];
    double sVal,sGrad0,sGrad1,sGrad2,sHess0,sHess1,sHess2,sHess4,sHess5,sHess8;
    
    complex_double *mod_coefs[64];
    intptr_t xs, ys, zs;
    int i, j, k, n, M;
    double *v,*g,*h, *p;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;
    
    tx = modf (ux, &ipartx);  int ix = (int) ipartx;
    ty = modf (uy, &iparty);  int iy = (int) iparty;
    tz = modf (uz, &ipartz);  int iz = (int) ipartz;

    a[0] = ( ( Ad[0]  * tx + Ad[1] ) * tx + Ad[2] ) * tx + Ad[3]; 
    a[1] = ( ( Ad[4]  * tx + Ad[5] ) * tx + Ad[6] ) * tx + Ad[7]; 
    a[2] = ( ( Ad[8]  * tx + Ad[9] ) * tx + Ad[10] ) * tx + Ad[11]; 
    a[3] = ( ( Ad[12] * tx + Ad[13] ) * tx + Ad[14] ) * tx + Ad[15]; 
    da[0] = ( ( dAd[0]  * tx + dAd[1] ) * tx + dAd[2] ) * tx + dAd[3]; 
    da[1] = ( ( dAd[4]  * tx + dAd[5] ) * tx + dAd[6] ) * tx + dAd[7]; 
    da[2] = ( ( dAd[8]  * tx + dAd[9] ) * tx + dAd[10] ) * tx + dAd[11]; 
    da[3] = ( ( dAd[12] * tx + dAd[13] ) * tx + dAd[14] ) * tx + dAd[15]; 
    d2a[0] = ( ( d2Ad[0]  * tx + d2Ad[1] ) * tx + d2Ad[2] ) * tx + d2Ad[3]; 
    d2a[1] = ( ( d2Ad[4]  * tx + d2Ad[5] ) * tx + d2Ad[6] ) * tx + d2Ad[7]; 
    d2a[2] = ( ( d2Ad[8]  * tx + d2Ad[9] ) * tx + d2Ad[10] ) * tx + d2Ad[11]; 
    d2a[3] = ( ( d2Ad[12] * tx + d2Ad[13] ) * tx + d2Ad[14] ) * tx + d2Ad[15]; 

    
    b[0] = ( ( Ad[0]  * ty + Ad[1] ) * ty + Ad[2] ) * ty + Ad[3]; 
    b[1] = ( ( Ad[4]  * ty + Ad[5] ) * ty + Ad[6] ) * ty + Ad[7]; 
    b[2] = ( ( Ad[8]  * ty + Ad[9] ) * ty + Ad[10] ) * ty + Ad[11]; 
    b[3] = ( ( Ad[12] * ty + Ad[13] ) * ty + Ad[14] ) * ty + Ad[15]; 
    db[0] = ( ( dAd[0]  * ty + dAd[1] ) * ty + dAd[2] ) * ty + dAd[3]; 
    db[1] = ( ( dAd[4]  * ty + dAd[5] ) * ty + dAd[6] ) * ty + dAd[7]; 
    db[2] = ( ( dAd[8]  * ty + dAd[9] ) * ty + dAd[10] ) * ty + dAd[11]; 
    db[3] = ( ( dAd[12] * ty + dAd[13] ) * ty + dAd[14] ) * ty + dAd[15]; 
    d2b[0] = ( ( d2Ad[0]  * ty + d2Ad[1] ) * ty + d2Ad[2] ) * ty + d2Ad[3]; 
    d2b[1] = ( ( d2Ad[4]  * ty + d2Ad[5] ) * ty + d2Ad[6] ) * ty + d2Ad[7]; 
    d2b[2] = ( ( d2Ad[8]  * ty + d2Ad[9] ) * ty + d2Ad[10] ) * ty + d2Ad[11]; 
    d2b[3] = ( ( d2Ad[12] * ty + d2Ad[13] ) * ty + d2Ad[14] ) * ty + d2Ad[15]; 
   
   
 
    c[0] = ( ( Ad[0]  * tz + Ad[1] ) * tz + Ad[2] ) * tz + Ad[3]; 
    c[1] = ( ( Ad[4]  * tz + Ad[5] ) * tz + Ad[6] ) * tz + Ad[7]; 
    c[2] = ( ( Ad[8]  * tz + Ad[9] ) * tz + Ad[10] ) * tz + Ad[11]; 
    c[3] = ( ( Ad[12] * tz + Ad[13] ) * tz + Ad[14] ) * tz + Ad[15]; 
    dc[0] = ( ( dAd[0]  * tz + dAd[1] ) * tz + dAd[2] ) * tz + dAd[3]; 
    dc[1] = ( ( dAd[4]  * tz + dAd[5] ) * tz + dAd[6] ) * tz + dAd[7]; 
    dc[2] = ( ( dAd[8]  * tz + dAd[9] ) * tz + dAd[10] ) * tz + dAd[11]; 
    dc[3] = ( ( dAd[12] * tz + dAd[13] ) * tz + dAd[14] ) * tz + dAd[15]; 
    d2c[0] = ( ( d2Ad[0]  * tz + d2Ad[1] ) * tz + d2Ad[2] ) * tz + d2Ad[3]; 
    d2c[1] = ( ( d2Ad[4]  * tz + d2Ad[5] ) * tz + d2Ad[6] ) * tz + d2Ad[7]; 
    d2c[2] = ( ( d2Ad[8]  * tz + d2Ad[9] ) * tz + d2Ad[10] ) * tz + d2Ad[11]; 
    d2c[3] = ( ( d2Ad[12] * tz + d2Ad[13] ) * tz + d2Ad[14] ) * tz + d2Ad[15]; 



    xs = spline->x_stride;
    ys = spline->y_stride;
    zs = spline->z_stride;

    double prefact[6];
    n = 0;
    for ( i=0; i<4; i++)
      for ( j=0; j<4; j++){
        prefact[0]=a[i]*b[j];
        prefact[1]=da[i]*b[j];
        prefact[2]=a[i]*db[j];
        prefact[3]=d2a[i]*b[j];
        prefact[4]=da[i]*db[j];
        prefact[5]=a[i]*d2b[j];
        for ( k=0; k<4; k++) {
	  dVal[n]   = prefact[0] * c[k];
	  dGrad0[n] = prefact[1] * c[k];
	  dGrad1[n] = prefact[2] * c[k];
	  dGrad2[n] = prefact[0] * dc[k];
	  dHess0[n] = prefact[3] * c[k];
	  dHess1[n] = prefact[4] * c[k];
	  dHess2[n] = prefact[1] * dc[k];
	  dHess4[n] = prefact[5] * c[k];
	  dHess5[n] = prefact[2] * dc[k];
	  dHess8[n] = prefact[0] * d2c[k];
          mod_coefs[n] = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	  n = n + 1;
        }
    }

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }

    M = 2 * ( spline -> num_splines );
    v = (double *)&vals[0];
    g = (double *)&grads[0];
    h = (double *)&hess[0];
    int rem = M % 4;
    for ( i = 0; i < 64; i++ )
    {  
       sVal = dVal[i];
       sGrad0=dGrad0[i];
       sGrad1=dGrad1[i];
       sGrad2=dGrad2[i];
       sHess0=dHess0[i];
       sHess1=dHess1[i];
       sHess2=dHess2[i];
       sHess4=dHess4[i];
       sHess5=dHess5[i];
       sHess8=dHess8[i];
       p = (double *)mod_coefs[i];

       for ( n = 0; n < M - rem ; n = n + 4 ) 
       {
          
          double aVal0, aVal1, aVal2, aVal3;
          double b0, b1, b2, b3;

          double aGrad000, aGrad001, aGrad002, aGrad003;
          double aGrad010, aGrad011, aGrad012, aGrad013;
          double aGrad020, aGrad021, aGrad022, aGrad023;

          double aHess000, aHess001, aHess002, aHess003;
          double aHess010, aHess011, aHess012, aHess013;
          double aHess020, aHess021, aHess022, aHess023;
          double aHess040, aHess041, aHess042, aHess043;
          double aHess050, aHess051, aHess052, aHess053;
          double aHess080, aHess081, aHess082, aHess083;

          b0 = p[n];
          b1 = p[n+1];
          b2 = p[n+2];
          b3 = p[n+3];


          //Val 
          aVal0 = v[n];
          aVal1 = v[n+1];
          aVal2 = v[n+2];
          aVal3 = v[n+3];

          aVal0 = aVal0 + sVal * b0;
          aVal1 = aVal1 + sVal * b1;
          aVal2 = aVal2 + sVal * b2;
          aVal3 = aVal3 + sVal * b3;

          v[n  ] = aVal0;
          v[n+1] = aVal1;
          v[n+2] = aVal2;
          v[n+3] = aVal3;

          //

          //Grad0
          aGrad000 = g[3*(n+0)+0];
          aGrad001 = g[3*(n+0)+1];
          aGrad010 = g[3*(n+0)+2];
          aGrad011 = g[3*(n+0)+3];

          aGrad000 = aGrad000 + sGrad0 * b0;
          aGrad001 = aGrad001 + sGrad0 * b1;
          aGrad010 = aGrad010 + sGrad1 * b0;
          aGrad011 = aGrad011 + sGrad1 * b1;

          g[3*(n+0)+0] = aGrad000;
          g[3*(n+0)+1] = aGrad001;
          g[3*(n+0)+2] = aGrad010;
          g[3*(n+0)+3] = aGrad011;

//        //Grad1
          aGrad020 = g[3*(n+0)+4];
          aGrad021 = g[3*(n+0)+5];
          aGrad002 = g[3*(n+2)+0];
          aGrad003 = g[3*(n+2)+1];
       

          aGrad020 = aGrad020 + sGrad2 * b0;
          aGrad021 = aGrad021 + sGrad2 * b1;
          aGrad002 = aGrad002 + sGrad0 * b2;
          aGrad003 = aGrad003 + sGrad0 * b3;

          g[3*(n+0)+4] = aGrad020;
          g[3*(n+0)+5] = aGrad021;
          g[3*(n+2)+0] = aGrad002;
          g[3*(n+2)+1] = aGrad003;

          //Grad2
          aGrad012 = g[3*(n+2)+2];
          aGrad013 = g[3*(n+2)+3];
          aGrad022 = g[3*(n+2)+4];
          aGrad023 = g[3*(n+2)+5];
          
          aGrad012 = aGrad012 + sGrad1 * b2;
          aGrad013 = aGrad013 + sGrad1 * b3;
          aGrad022 = aGrad022 + sGrad2 * b2;
          aGrad023 = aGrad023 + sGrad2 * b3;
 
          g[3*(n+2)+2] = aGrad012;
          g[3*(n+2)+3] = aGrad013;
          g[3*(n+2)+4] = aGrad022;
          g[3*(n+2)+5] = aGrad023;



          //hess0
          aHess000 = h[9*(n+0)+0];
          aHess001 = h[9*(n+0)+1];
          aHess010 = h[9*(n+0)+2];
          aHess011 = h[9*(n+0)+3];

          aHess000 = aHess000 + sHess0 * b0;
          aHess001 = aHess001 + sHess0 * b1;
          aHess010 = aHess010 + sHess1 * b0;
          aHess011 = aHess011 + sHess1 * b1;

 
          h[9*(n+0)+0] = aHess000;
          h[9*(n+0)+1] = aHess001;
          h[9*(n+0)+2] = aHess010;
          h[9*(n+0)+3] = aHess011;

          
          
          //hess1
          aHess020 = h[9*(n+0)+4];
          aHess021 = h[9*(n+0)+5];
          aHess020 = aHess020 + sHess2 * b0;
          aHess021 = aHess021 + sHess2 * b1;
          h[9*(n+0)+4] = aHess020;
          h[9*(n+0)+5] = aHess021;



        // Hess2 
          aHess040 = h[9*(n+0)+8];
          aHess041 = h[9*(n+0)+9];
          aHess050 = h[9*(n+0)+10];
          aHess051 = h[9*(n+0)+11];


          aHess040 = aHess040 + sHess4 * b0;
          aHess041 = aHess041 + sHess4 * b1;
          aHess050 = aHess050 + sHess5 * b0;
          aHess051 = aHess051 + sHess5 * b1;

          h[9*(n+0)+8] = aHess040;
          h[9*(n+0)+9] = aHess041;
          h[9*(n+0)+10] = aHess050;
          h[9*(n+0)+11] = aHess051;


          //hess3

          aHess080 = h[9*(n+0)+16];
          aHess081 = h[9*(n+0)+17];
          aHess002 = h[9*(n+2)+0];
          aHess003 = h[9*(n+2)+1];

          aHess080 = aHess080 + sHess8 * b0;
          aHess081 = aHess081 + sHess8 * b1;
          aHess002 = aHess002 + sHess0 * b2;
          aHess003 = aHess003 + sHess0 * b3;


          h[9*(n+0)+16] = aHess080;
          h[9*(n+0)+17] = aHess081;
          h[9*(n+2)+0] = aHess002;
          h[9*(n+2)+1] = aHess003;


         //Hess4 
          aHess012 = h[9*(n+2)+2];
          aHess013 = h[9*(n+2)+3];
          aHess022 = h[9*(n+2)+4];
          aHess023 = h[9*(n+2)+5];

          aHess012 = aHess012 + sHess1 * b2;
          aHess013 = aHess013 + sHess1 * b3;
          aHess022 = aHess022 + sHess2 * b2;
          aHess023 = aHess023 + sHess2 * b3;

          h[9*(n+2)+2] = aHess012;
          h[9*(n+2)+3] = aHess013;
          h[9*(n+2)+4] = aHess022;
          h[9*(n+2)+5] = aHess023;

          //Hess5

          aHess042 = h[9*(n+2)+8];
          aHess043 = h[9*(n+2)+9];
          aHess052 = h[9*(n+2)+10];
          aHess053 = h[9*(n+2)+11];
          

          aHess042 = aHess042 + sHess4 * b2;
          aHess043 = aHess043 + sHess4 * b3;
          aHess052 = aHess052 + sHess5 * b2;
          aHess053 = aHess053 + sHess5 * b3;
 

          h[9*(n+2)+8] = aHess042;
          h[9*(n+2)+9] = aHess043;
          h[9*(n+2)+10] = aHess052;
          h[9*(n+2)+11] = aHess053;

          //Hess6
          
          aHess082 = h[9*(n+2)+16];
          aHess083 = h[9*(n+2)+17];
          aHess082 = aHess082 + sHess8 * b2;
          aHess083 = aHess083 + sHess8 * b3;
          h[9*(n+2)+16] = aHess082;
          h[9*(n+2)+17] = aHess083;
       }

       for ( k = n; k < M; k++){
           v[k] = v[k] + sVal * p[k];
       }
          for ( k = n; k < M; k=k+2 ){
           g[3*k+0] = g[3*k+0] + sGrad0 * p[k];
           g[3*k+1] = g[3*k+1] + sGrad0 * p[k+1];
           g[3*k+2] = g[3*k+2] + sGrad1 * p[k];
           g[3*k+3] = g[3*k+3] + sGrad1 * p[k+1];
           g[3*k+4] = g[3*k+4] + sGrad2 * p[k];
           g[3*k+5] = g[3*k+5] + sGrad2 * p[k+1];
           h[9*k+0] = h[9*k+0] + sHess0 * p[k];
           h[9*k+1] = h[9*k+1] + sHess0 * p[k+1];
           h[9*k+2] = h[9*k+2] + sHess1 * p[k];
           h[9*k+3] = h[9*k+3] + sHess1 * p[k+1];
           h[9*k+4] = h[9*k+4] + sHess2 * p[k];
           h[9*k+5] = h[9*k+5] + sHess2 * p[k+1];
           h[9*k+8] = h[9*k+8] + sHess4 * p[k];
           h[9*k+9] = h[9*k+9] + sHess4 * p[k+1];
           h[9*k+10] = h[9*k+10] + sHess5 * p[k];
           h[9*k+11] = h[9*k+11] + sHess5 * p[k+1];
           h[9*k+16] = h[9*k+16] + sHess8 * p[k];
           h[9*k+17] = h[9*k+17] + sHess8 * p[k+1];
}
    } 

    dxInv = spline->x_grid.delta_inv;
    dyInv = spline->y_grid.delta_inv;
    dzInv = spline->z_grid.delta_inv; 
 
    dxInvdxInv = dxInv*dxInv;
    dyInvdyInv = dyInv*dyInv;
    dzInvdzInv = dzInv*dzInv;
    dxInvdyInv = dxInv*dyInv;
    dxInvdzInv = dxInv*dzInv;
    dyInvdzInv = dyInv*dzInv;
 
    for (int n=0; n<spline->num_splines; n++) {
        grads[3*n+0] *= dxInv;
        grads[3*n+1] *= dyInv;
        grads[3*n+2] *= dzInv;
        hess[9*n+0] *= dxInvdxInv;
        hess[9*n+4] *= dyInvdyInv;
        hess[9*n+8] *= dzInvdzInv;
        hess[9*n+1] *= dxInvdyInv;
        hess[9*n+2] *= dxInvdzInv;
        hess[9*n+5] *= dyInvdzInv;
        // Copy hessian elements into lower half of 3x3 matrix
        hess[9*n+3] = hess[9*n+1];
        hess[9*n+6] = hess[9*n+2];
        hess[9*n+7] = hess[9*n+5];

     }
}
void
eval_multi_UBspline_3d_z_vg (const multi_UBspline_3d_z *spline,
			      double x, double y, double z,
			      complex_double* restrict vals,
			     complex_double* restrict grads)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 0]*tpz[0] + dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 4]*tpz[0] + dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 8]*tpz[0] + dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[12]*tpz[0] + dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) 
      for (int k=0; k<4; k++) {
	double abc = a[i]*b[j]*c[k];
	double dabc[3];
	dabc[0] = da[i]* b[j]* c[k];
	dabc[1] =  a[i]*db[j]* c[k];
	dabc[2] =  a[i]* b[j]*dc[k];

	complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	for (int n=0; n<spline->num_splines; n++) {
	  vals[n]      +=   abc   *coefs[n];
	  grads[3*n+0] +=  dabc[0]*coefs[n];
	  grads[3*n+1] +=  dabc[1]*coefs[n];
	  grads[3*n+2] +=  dabc[2]*coefs[n];
	}
      }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv; 
  for (int n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
  }
}



void
eval_multi_UBspline_3d_z_vgl (const multi_UBspline_3d_z *spline,
			      double x, double y, double z,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict lapl)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

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

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 0]*tpz[0] + dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 4]*tpz[0] + dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 8]*tpz[0] + dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[12]*tpz[0] + dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 0]*tpz[0] + d2Ad[ 1]*tpz[1] + d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 4]*tpz[0] + d2Ad[ 5]*tpz[1] + d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[ 8]*tpz[0] + d2Ad[ 9]*tpz[1] + d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[12]*tpz[0] + d2Ad[13]*tpz[1] + d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  //complex_double lapl3[3*spline->num_splines];
  complex_double* restrict lapl3 = spline->lapl3;
  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    lapl3[3*n+0] = lapl3[3*n+1] = lapl3[3*n+2] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) 
      for (int k=0; k<4; k++) {
	double abc = a[i]*b[j]*c[k];
	double dabc[3], d2abc[3];
	dabc[0] = da[i]* b[j]* c[k];
	dabc[1] =  a[i]*db[j]* c[k];
	dabc[2] =  a[i]* b[j]*dc[k];
	d2abc[0] = d2a[i]*  b[j]*  c[k];
	d2abc[1] =   a[i]*d2b[j]*  c[k];
	d2abc[2] =   a[i]*  b[j]*d2c[k];

	complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	for (int n=0; n<spline->num_splines; n++) {
	  vals[n]      +=   abc   *coefs[n];
	  grads[3*n+0] +=  dabc[0]*coefs[n];
	  grads[3*n+1] +=  dabc[1]*coefs[n];
	  grads[3*n+2] +=  dabc[2]*coefs[n];
	  lapl3[3*n+0] += d2abc[0]*coefs[n];
	  lapl3[3*n+1] += d2abc[1]*coefs[n];
	  lapl3[3*n+2] += d2abc[2]*coefs[n];
	}
      }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv; 
  for (int n=0; n<spline->num_splines; n++) {
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
eval_multi_UBspline_3d_z_vghgh (const multi_UBspline_3d_z *spline,
               double x, double y, double z,
               complex_double* restrict vals,
               complex_double* restrict grads,
               complex_double* restrict hess,
               complex_double* restrict gradhess)    
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4],
    d3a[4], d3b[4], d3c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  complex_double* restrict coefs = spline->coefs;

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
  d3a[0] = (/*d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] +*/ d3Ad[ 3]*tpx[3]);
  d3a[1] = (/*d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] +*/ d3Ad[ 7]*tpx[3]);
  d3a[2] = (/*d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] +*/ d3Ad[11]*tpx[3]);
  d3a[3] = (/*d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] +*/ d3Ad[15]*tpx[3]);
  
  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);
  d3b[0] = (/*d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] +*/ d3Ad[ 3]*tpy[3]);
  d3b[1] = (/*d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] +*/ d3Ad[ 7]*tpy[3]);
  d3b[2] = (/*d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] +*/ d3Ad[11]*tpy[3]);
  d3b[3] = (/*d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] +*/ d3Ad[15]*tpy[3]);

  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 0]*tpz[0] + dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 4]*tpz[0] + dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 8]*tpz[0] + dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[12]*tpz[0] + dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 0]*tpz[0] + d2Ad[ 1]*tpz[1] + d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 4]*tpz[0] + d2Ad[ 5]*tpz[1] + d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[ 8]*tpz[0] + d2Ad[ 9]*tpz[1] + d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[12]*tpz[0] + d2Ad[13]*tpz[1] + d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);
  d3c[0] = (/*d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] +*/ d3Ad[ 3]*tpz[3]);
  d3c[1] = (/*d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] +*/ d3Ad[ 7]*tpz[3]);
  d3c[2] = (/*d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] +*/ d3Ad[11]*tpz[3]);
  d3c[3] = (/*d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] +*/ d3Ad[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
    for (int i=0; i<27; i++)
      gradhess[27*n+i] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) 
      for (int k=0; k<4; k++) {
   double abc = a[i]*b[j]*c[k];
   double dabc[3], d2abc[6], d3abc[10];
   dabc[0] = da[i]* b[j]* c[k];
   dabc[1] =  a[i]*db[j]* c[k];
   dabc[2] =  a[i]* b[j]*dc[k];
   d2abc[0] = d2a[i]*  b[j]*  c[k];
   d2abc[1] =  da[i]* db[j]*  c[k];
   d2abc[2] =  da[i]*  b[j]* dc[k];
   d2abc[3] =   a[i]*d2b[j]*  c[k];
   d2abc[4] =   a[i]* db[j]* dc[k];
   d2abc[5] =   a[i]*  b[j]*d2c[k];
   
   d3abc[0] = d3a[i]*  b[j]*  c[k];
   d3abc[1] = d2a[i]* db[j]*  c[k];
   d3abc[2] = d2a[i]*  b[j]* dc[k];
   d3abc[3] =  da[i]*d2b[j]*  c[k];
   d3abc[4] =  da[i]* db[j]* dc[k];
   d3abc[5] =  da[i]*  b[j]*d2c[k];   
   d3abc[6] =   a[i]*d3b[j]*  c[k];
   d3abc[7] =   a[i]*d2b[j]* dc[k];
   d3abc[8] =   a[i]* db[j]*d2c[k];     
   d3abc[9] =   a[i]*  b[j]*d3c[k];

   complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
   for (int n=0; n<spline->num_splines; n++) {
     vals[n]      +=   abc   *coefs[n];
     grads[3*n+0] +=  dabc[0]*coefs[n];
     grads[3*n+1] +=  dabc[1]*coefs[n];
     grads[3*n+2] +=  dabc[2]*coefs[n];
     hess [9*n+0] += d2abc[0]*coefs[n];
     hess [9*n+1] += d2abc[1]*coefs[n];
     hess [9*n+2] += d2abc[2]*coefs[n];
     hess [9*n+4] += d2abc[3]*coefs[n];
     hess [9*n+5] += d2abc[4]*coefs[n];
     hess [9*n+8] += d2abc[5]*coefs[n];
     
     gradhess [27*n+0 ] += d3abc[0]*coefs[n];
     gradhess [27*n+1 ] += d3abc[1]*coefs[n];
     gradhess [27*n+2 ] += d3abc[2]*coefs[n];
     gradhess [27*n+4 ] += d3abc[3]*coefs[n];
     gradhess [27*n+5 ] += d3abc[4]*coefs[n];
     gradhess [27*n+8 ] += d3abc[5]*coefs[n];
     gradhess [27*n+13] += d3abc[6]*coefs[n];
     gradhess [27*n+14] += d3abc[7]*coefs[n];
     gradhess [27*n+17] += d3abc[8]*coefs[n];
     gradhess [27*n+26] += d3abc[9]*coefs[n];
   }
      }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv; 
  for (int n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess [9*n+0] *= dxInv*dxInv;
    hess [9*n+4] *= dyInv*dyInv;
    hess [9*n+8] *= dzInv*dzInv;
    hess [9*n+1] *= dxInv*dyInv;
    hess [9*n+2] *= dxInv*dzInv;
    hess [9*n+5] *= dyInv*dzInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess [9*n+3] = hess[9*n+1];
    hess [9*n+6] = hess[9*n+2];
    hess [9*n+7] = hess[9*n+5];
    
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

#endif
