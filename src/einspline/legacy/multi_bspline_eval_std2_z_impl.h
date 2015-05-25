////////////////////////////////////////////////////////////////////////////
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

#ifdef BGQPX 
#include <builtins.h>
void
eval_multi_UBspline_3d_z (const multi_UBspline_3d_z *spline,
                          double x, double y, double z,
                          complex_double* restrict vals)
{
    double ux, uy, uz, prefactor, ipartx, iparty, ipartz, tx, ty, tz, a[4], b[4],c[4], d[64], s;
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
    int rem=M % 8; 
    for ( n = 0; n < M; n++ ) v[n] = 0.0;

    for ( i = 0; i < 64; i++ )
    {   
       s = d[i];
       __dcbt(&mod_coefs[i+3]);
       p = (double *)mod_coefs[i];
       vector4double t = { s, s, s, s };

       for ( n = 0, j = 0; n < M-rem; n = n + 8, j = j + 64 ) 
       {   


          vector4double f0, f1; 
//          vector4double g0, g1; 

          __dcbt(&p[n+32]);
          __dcbt(&v[n+8]);

          //g0 = vec_ld( j,    p );
          //g1 = vec_ld( j+32, p );

          vector4double  g0= {p[n+0],p[n+1],p[n+2],p[n+3]};
          vector4double  g1= {p[n+4],p[n+5],p[n+6],p[n+7]};

          f0 = vec_ld( j,    v );
          f1 = vec_ld( j+32, v );



          f0 = vec_madd( t, g0, f0 );
          f1 = vec_madd( t, g1, f1 );

          vec_st( f0, j,    v );
          vec_st( f1, j+32, v );

       }

       for ( k = n; k < M; k++ )
           v[k] = v[k] + s * p[k];

    }
}

void
eval_multi_UBspline_3d_z_vgh (const multi_UBspline_3d_z *spline,
			      double x, double y, double z,
			      complex_double* restrict vals,
			      complex_double* restrict grads,
			      complex_double* restrict hess)
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

  double ap[4],bp[4],cp[4];
  double dap[4],dbp[4],dcp[4];
  double d2ap[4],d2bp[4],d2cp[4];
  ap[0] = ( ( Ad[0]  * tx + Ad[1] ) * tx + Ad[2] ) * tx + Ad[3];
  ap[1] = ( ( Ad[4]  * tx + Ad[5] ) * tx + Ad[6] ) * tx + Ad[7];
  ap[2] = ( ( Ad[8]  * tx + Ad[9] ) * tx + Ad[10] ) * tx + Ad[11];
  ap[3] = ( ( Ad[12] * tx + Ad[13] ) * tx + Ad[14] ) * tx + Ad[15]; 
  dap[0] = ( ( dAd[0]  * tx + dAd[1] ) * tx + dAd[2] ) * tx + dAd[3]; 
  dap[1] = ( ( dAd[4]  * tx + dAd[5] ) * tx + dAd[6] ) * tx + dAd[7]; 
  dap[2] = ( ( dAd[8]  * tx + dAd[9] ) * tx + dAd[10] ) * tx + dAd[11]; 
  dap[3] = ( ( dAd[12] * tx + dAd[13] ) * tx + dAd[14] ) * tx + dAd[15]; 
  d2ap[0] = ( ( d2Ad[0]  * tx + d2Ad[1] ) * tx + d2Ad[2] ) * tx + d2Ad[3]; 
  d2ap[1] = ( ( d2Ad[4]  * tx + d2Ad[5] ) * tx + d2Ad[6] ) * tx + d2Ad[7]; 
  d2ap[2] = ( ( d2Ad[8]  * tx + d2Ad[9] ) * tx + d2Ad[10] ) * tx + d2Ad[11]; 
  d2ap[3] = ( ( d2Ad[12] * tx + d2Ad[13] ) * tx + d2Ad[14] ) * tx + d2Ad[15]; 
    
  bp[0] = ( ( Ad[0]  * ty + Ad[1] ) * ty + Ad[2] ) * ty + Ad[3]; 
  bp[1] = ( ( Ad[4]  * ty + Ad[5] ) * ty + Ad[6] ) * ty + Ad[7]; 
  bp[2] = ( ( Ad[8]  * ty + Ad[9] ) * ty + Ad[10] ) * ty + Ad[11]; 
  bp[3] = ( ( Ad[12] * ty + Ad[13] ) * ty + Ad[14] ) * ty + Ad[15]; 
  dbp[0] = ( ( dAd[0]  * ty + dAd[1] ) * ty + dAd[2] ) * ty + dAd[3]; 
  dbp[1] = ( ( dAd[4]  * ty + dAd[5] ) * ty + dAd[6] ) * ty + dAd[7]; 
  dbp[2] = ( ( dAd[8]  * ty + dAd[9] ) * ty + dAd[10] ) * ty + dAd[11]; 
  dbp[3] = ( ( dAd[12] * ty + dAd[13] ) * ty + dAd[14] ) * ty + dAd[15]; 
  d2bp[0] = ( ( d2Ad[0]  * ty + d2Ad[1] ) * ty + d2Ad[2] ) * ty + d2Ad[3]; 
  d2bp[1] = ( ( d2Ad[4]  * ty + d2Ad[5] ) * ty + d2Ad[6] ) * ty + d2Ad[7]; 
  d2bp[2] = ( ( d2Ad[8]  * ty + d2Ad[9] ) * ty + d2Ad[10] ) * ty + d2Ad[11]; 
  d2bp[3] = ( ( d2Ad[12] * ty + d2Ad[13] ) * ty + d2Ad[14] ) * ty + d2Ad[15]; 
    
  cp[0] = ( ( Ad[0]  * tz + Ad[1] ) * tz + Ad[2] ) * tz + Ad[3]; 
  cp[1] = ( ( Ad[4]  * tz + Ad[5] ) * tz + Ad[6] ) * tz + Ad[7]; 
  cp[2] = ( ( Ad[8]  * tz + Ad[9] ) * tz + Ad[10] ) * tz + Ad[11]; 
  cp[3] = ( ( Ad[12] * tz + Ad[13] ) * tz + Ad[14] ) * tz + Ad[15]; 
  dcp[0] = ( ( dAd[0]  * tz + dAd[1] ) * tz + dAd[2] ) * tz + dAd[3]; 
  dcp[1] = ( ( dAd[4]  * tz + dAd[5] ) * tz + dAd[6] ) * tz + dAd[7]; 
  dcp[2] = ( ( dAd[8]  * tz + dAd[9] ) * tz + dAd[10] ) * tz + dAd[11]; 
  dcp[3] = ( ( dAd[12] * tz + dAd[13] ) * tz + dAd[14] ) * tz + dAd[15]; 
  d2cp[0] = ( ( d2Ad[0]  * tz + d2Ad[1] ) * tz + d2Ad[2] ) * tz + d2Ad[3]; 
  d2cp[1] = ( ( d2Ad[4]  * tz + d2Ad[5] ) * tz + d2Ad[6] ) * tz + d2Ad[7]; 
  d2cp[2] = ( ( d2Ad[8]  * tz + d2Ad[9] ) * tz + d2Ad[10] ) * tz + d2Ad[11]; 
  d2cp[3] = ( ( d2Ad[12] * tz + d2Ad[13] ) * tz + d2Ad[14] ) * tz + d2Ad[15]; 



  vector4double a= vec_ld(sizeof(double)*0,ap);
  vector4double b= vec_ld(sizeof(double)*0,bp);
  vector4double c= vec_ld(sizeof(double)*0,cp);
  vector4double da= vec_ld(sizeof(double)*0,dap);
  vector4double db= vec_ld(sizeof(double)*0,dbp);
  vector4double dc= vec_ld(sizeof(double)*0,dcp);
  vector4double d2a= vec_ld(sizeof(double)*0,d2ap);
  vector4double d2b= vec_ld(sizeof(double)*0,d2bp);
  vector4double d2c= vec_ld(sizeof(double)*0,d2cp);


  vector4double mantiss_1 = vec_gpci(0011);
  vector4double mantiss_2 = vec_gpci(02233);

  vector4double splatc_0_1=vec_perm(c,c,mantiss_1);
  vector4double splatc_2_3=vec_perm(c,c,mantiss_2);

  vector4double splatdc_0_1=vec_perm(dc,dc,mantiss_1);
  vector4double splatdc_2_3=vec_perm(dc,dc,mantiss_2);

  vector4double splatd2c_0_1=vec_perm(d2c,d2c,mantiss_1);
  vector4double splatd2c_2_3=vec_perm(d2c,d2c,mantiss_2);


  vector4double b_0=vec_splat(b,0);
  vector4double b_1=vec_splat(b,1);
  vector4double b_2=vec_splat(b,2);
  vector4double b_3=vec_splat(b,3);
   
  vector4double db_0=vec_splat(db,0);
  vector4double db_1=vec_splat(db,1);
  vector4double db_2=vec_splat(db,2);
  vector4double db_3=vec_splat(db,3);
   
  vector4double d2b_0=vec_splat(d2b,0);
  vector4double d2b_1=vec_splat(d2b,1);
  vector4double d2b_2=vec_splat(d2b,2);
  vector4double d2b_3=vec_splat(d2b,3);
   

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;


  vector4double val_temp,grad0_temp, grad1_temp, grad2_temp, hess0_temp, hess1_temp, hess2_temp, hess4_temp, hess5_temp, hess8_temp; 


  complex_double* restrict coefs = spline->coefs + ix*xs + iy*ys + iz*zs;
  int num_splines = spline->num_splines;
  for (int n=0; n<num_splines; n++, coefs++) {

      vector4double val = {0.0};
      vector4double grad0 = {0.0};
      vector4double grad1 = {0.0};
      vector4double grad2 = {0.0};
      vector4double hess0 = {0.0};
      vector4double hess1 = {0.0};
      vector4double hess2 = {0.0};
      vector4double hess4 = {0.0};
      vector4double hess5 = {0.0};
      vector4double hess8 = {0.0};
      for (int i=0; i<4; i++) {
         //  First Coefficient 
         vector4double Coef0_0_0=vec_ld2(0,(double*)(coefs+(i*xs)));
         vector4double Coef0_0_1=vec_ld2(0,(double*)(coefs+(i*xs+zs)));
         vector4double Coef0_1_0=vec_sldw(Coef0_0_0,Coef0_0_1,2);

         vector4double Coef0_0_2=vec_ld2(0,(double*)(coefs+(i*xs+2*zs)));
         vector4double Coef0_0_3=vec_ld2(0,(double*)(coefs+(i*xs+3*zs)));
         vector4double Coef0_2_0=vec_sldw(Coef0_0_2,Coef0_0_3,2);


         vector4double Coef_Prime_0_0 =vec_mul(Coef0_1_0,splatc_0_1);
         vector4double Coef_Prime_0_1 =vec_mul(Coef0_2_0,splatc_2_3);
         vector4double Coef_Primedc_0_0 =vec_mul(Coef0_1_0,splatdc_0_1);
         vector4double Coef_Primedc_0_1 =vec_mul(Coef0_2_0,splatdc_2_3);
         vector4double Coef_Primed2c_0_0 =vec_mul(Coef0_1_0,splatd2c_0_1);
         vector4double Coef_Primed2c_0_1 =vec_mul(Coef0_2_0,splatd2c_2_3);


         vector4double Coef0_0= vec_add(Coef_Prime_0_0,Coef_Prime_0_1);
         vector4double Coefdc0_0= vec_add(Coef_Primedc_0_0,Coef_Primedc_0_1);
         vector4double Coefd2c0_0= vec_add(Coef_Primed2c_0_0,Coef_Primed2c_0_1);

         vector4double Coef0=vec_mul(Coef0_0,b_0);
         vector4double Coefdb0=vec_mul(Coef0_0,db_0);
         vector4double Coefd2b0=vec_mul(Coef0_0,d2b_0);
         vector4double Coefdc0=vec_mul(Coefdc0_0,b_0);
         vector4double Coefdbdc0=vec_mul(Coefdc0_0,db_0);
         vector4double Coefd2c0=vec_mul(Coefd2c0_0,b_0);

         //  Second Coefficient 
         vector4double Coef1_0_0=vec_ld2(0,(double*)(coefs+(i*xs+ys)));
         vector4double Coef1_0_1=vec_ld2(0,(double*)(coefs+(i*xs+ys+zs)));
         vector4double Coef1_1_0=vec_sldw(Coef1_0_0,Coef1_0_1,2);

         vector4double Coef1_0_2=vec_ld2(0,(double*)(coefs+(i*xs+ys+2*zs)));
         vector4double Coef1_0_3=vec_ld2(0,(double*)(coefs+(i*xs+ys+3*zs)));
         vector4double Coef1_2_0=vec_sldw(Coef1_0_2,Coef1_0_3,2);


         vector4double Coef_Prime_1_0 =vec_mul(Coef1_1_0,splatc_0_1);
         vector4double Coef_Prime_1_1 =vec_mul(Coef1_2_0,splatc_2_3);
         vector4double Coef_Primedc_1_0 =vec_mul(Coef1_1_0,splatdc_0_1);
         vector4double Coef_Primedc_1_1 =vec_mul(Coef1_2_0,splatdc_2_3);
         vector4double Coef_Primed2c_1_0 =vec_mul(Coef1_1_0,splatd2c_0_1);
         vector4double Coef_Primed2c_1_1 =vec_mul(Coef1_2_0,splatd2c_2_3);


         vector4double Coef1_0= vec_add(Coef_Prime_1_0,Coef_Prime_1_1);
         vector4double Coefdc1_0= vec_add(Coef_Primedc_1_0,Coef_Primedc_1_1);
         vector4double Coefd2c1_0= vec_add(Coef_Primed2c_1_0,Coef_Primed2c_1_1);

         vector4double Coef1=vec_mul(Coef1_0,b_1);
         vector4double Coefdb1=vec_mul(Coef1_0,db_1);
         vector4double Coefd2b1=vec_mul(Coef1_0,d2b_1);
         vector4double Coefdc1=vec_mul(Coefdc1_0,b_1);
         vector4double Coefdbdc1=vec_mul(Coefdc1_0,db_1);
         vector4double Coefd2c1=vec_mul(Coefd2c1_0,b_1);

         //  Third Coefficient 
         vector4double Coef2_0_0=vec_ld2(0,(double*)(coefs+(i*xs+2*ys)));
         vector4double Coef2_0_1=vec_ld2(0,(double*)(coefs+(i*xs+2*ys+zs)));
         vector4double Coef2_1_0=vec_sldw(Coef2_0_0,Coef2_0_1,2);

         vector4double Coef2_0_2=vec_ld2(0,(double*)(coefs+(i*xs+2*ys+2*zs)));
         vector4double Coef2_0_3=vec_ld2(0,(double*)(coefs+(i*xs+2*ys+3*zs)));
         vector4double Coef2_2_0=vec_sldw(Coef2_0_2,Coef2_0_3,2);


         vector4double Coef_Prime_2_0 =vec_mul(Coef2_1_0,splatc_0_1);
         vector4double Coef_Prime_2_1 =vec_mul(Coef2_2_0,splatc_2_3);
         vector4double Coef_Primedc_2_0 =vec_mul(Coef2_1_0,splatdc_0_1);
         vector4double Coef_Primedc_2_1 =vec_mul(Coef2_2_0,splatdc_2_3);
         vector4double Coef_Primed2c_2_0 =vec_mul(Coef2_1_0,splatd2c_0_1);
         vector4double Coef_Primed2c_2_1 =vec_mul(Coef2_2_0,splatd2c_2_3);

         vector4double Coef2_0= vec_add(Coef_Prime_2_0,Coef_Prime_2_1);
         vector4double Coefdc2_0= vec_add(Coef_Primedc_2_0,Coef_Primedc_2_1);
         vector4double Coefd2c2_0= vec_add(Coef_Primed2c_2_0,Coef_Primed2c_2_1);

         vector4double Coef2=vec_mul(Coef2_0,b_2);
         vector4double Coefdb2=vec_mul(Coef2_0,db_2);
         vector4double Coefd2b2=vec_mul(Coef2_0,d2b_2);
         vector4double Coefdc2=vec_mul(Coefdc2_0,b_2);
         vector4double Coefdbdc2=vec_mul(Coefdc2_0,db_2);
         vector4double Coefd2c2=vec_mul(Coefd2c2_0,b_2);


         //  Fourth Coefficient 
         vector4double Coef3_0_0=vec_ld2(0,(double*)(coefs+(i*xs+3*ys)));
         vector4double Coef3_0_1=vec_ld2(0,(double*)(coefs+(i*xs+3*ys+zs)));
         vector4double Coef3_1_0=vec_sldw(Coef3_0_0,Coef3_0_1,2);

         vector4double Coef3_0_2=vec_ld2(0,(double*)(coefs+(i*xs+3*ys+2*zs)));
         vector4double Coef3_0_3=vec_ld2(0,(double*)(coefs+(i*xs+3*ys+3*zs)));
         vector4double Coef3_2_0=vec_sldw(Coef3_0_2,Coef3_0_3,2);


         vector4double Coef_Prime_3_0 =vec_mul(Coef3_1_0,splatc_0_1);
         vector4double Coef_Prime_3_1 =vec_mul(Coef3_2_0,splatc_2_3);
         vector4double Coef_Primedc_3_0 =vec_mul(Coef3_1_0,splatdc_0_1);
         vector4double Coef_Primedc_3_1 =vec_mul(Coef3_2_0,splatdc_2_3);
         vector4double Coef_Primed2c_3_0 =vec_mul(Coef3_1_0,splatd2c_0_1);
         vector4double Coef_Primed2c_3_1 =vec_mul(Coef3_2_0,splatd2c_2_3);


         vector4double Coef3_0= vec_add(Coef_Prime_3_0,Coef_Prime_3_1);
         vector4double Coefdc3_0= vec_add(Coef_Primedc_3_0,Coef_Primedc_3_1);
         vector4double Coefd2c3_0= vec_add(Coef_Primed2c_3_0,Coef_Primed2c_3_1);

         vector4double Coef3=vec_mul(Coef3_0,b_3);
         vector4double Coefdb3=vec_mul(Coef3_0,db_3);
         vector4double Coefd2b3=vec_mul(Coef3_0,d2b_3);
         vector4double Coefdc3=vec_mul(Coefdc3_0,b_3);
         vector4double Coefdbdc3=vec_mul(Coefdc3_0,db_3);
         vector4double Coefd2c3=vec_mul(Coefd2c3_0,b_3);

         //Summing all coefficients
         vector4double Coef_sum0=vec_add(Coef0,Coef1);
         vector4double Coef_sum1=vec_add(Coef2,Coef3);
         vector4double Coef=vec_add(Coef_sum0,Coef_sum1);

         vector4double Coefdb_sum0=vec_add(Coefdb0,Coefdb1);
         vector4double Coefdb_sum1=vec_add(Coefdb2,Coefdb3);
         vector4double Coefdb=vec_add(Coefdb_sum0,Coefdb_sum1);

         vector4double Coefd2b_sum0=vec_add(Coefd2b0,Coefd2b1);
         vector4double Coefd2b_sum1=vec_add(Coefd2b2,Coefd2b3);
         vector4double Coefd2b=vec_add(Coefd2b_sum0,Coefd2b_sum1);

         vector4double Coefdc_sum0=vec_add(Coefdc0,Coefdc1);
         vector4double Coefdc_sum1=vec_add(Coefdc2,Coefdc3);
         vector4double Coefdc=vec_add(Coefdc_sum0,Coefdc_sum1);

         vector4double Coefd2c_sum0=vec_add(Coefd2c0,Coefd2c1);
         vector4double Coefd2c_sum1=vec_add(Coefd2c2,Coefd2c3);
         vector4double Coefd2c=vec_add(Coefd2c_sum0,Coefd2c_sum1);


         vector4double Coefdbdc_sum0=vec_add(Coefdbdc0,Coefdbdc1);
         vector4double Coefdbdc_sum1=vec_add(Coefdbdc2,Coefdbdc3);
         vector4double Coefdbdc=vec_add(Coefdbdc_sum0,Coefdbdc_sum1);



         // Multiplying Coeffs by A[i] 
         vector4double asplat = vec_splats(a[i]);
         vector4double dasplat = vec_splats(da[i]);
         vector4double d2asplat = vec_splats(d2a[i]);
         
         val_temp=vec_mul(Coef,asplat);
         grad0_temp=vec_mul(Coef,dasplat);
         grad1_temp=vec_mul(Coefdb,asplat);
         grad2_temp=vec_mul(Coefdc,asplat);
    
         hess0_temp=vec_mul(Coef,d2asplat);
         hess1_temp=vec_mul(Coefdb,dasplat);
         hess2_temp=vec_mul(Coefdc,dasplat);
         hess4_temp=vec_mul(Coefd2b,asplat);
         hess5_temp=vec_mul(Coefdbdc,asplat);  
         hess8_temp=vec_mul(Coefd2c,asplat);

         //Sum of coefficients over loop;
         val=vec_add(val_temp,val);

         grad0=vec_add(grad0_temp,grad0);
         grad1=vec_add(grad1_temp,grad1);
         grad2=vec_add(grad2_temp,grad2);

         hess0=vec_add(hess0_temp,hess0);
         hess1=vec_add(hess1_temp,hess1);
         hess2=vec_add(hess2_temp,hess2);
         hess4=vec_add(hess4_temp,hess4);
         hess5=vec_add(hess5_temp,hess5);
         hess8=vec_add(hess8_temp,hess8);

      }
      complex_double temp_val = (val[0]+val[2])+(val[1]+val[3])*_Complex_I;
      vals[n] = temp_val; 

      complex_double temp_grads0 = (grad0[0]+grad0[2])+(grad0[1]+grad0[3])*_Complex_I;
      grads[3*n+0] = temp_grads0; 
      complex_double temp_grads1 = (grad1[0]+grad1[2])+(grad1[1]+grad1[3])*_Complex_I;
      grads[3*n+1] = temp_grads1; 
      complex_double temp_grads2 = (grad2[0]+grad2[2])+(grad2[1]+grad2[3])*_Complex_I;
      grads[3*n+2] = temp_grads2; 

      complex_double temp_hess0 = (hess0[0]+hess0[2])+(hess0[1]+hess0[3])*_Complex_I;
      hess[9*n+0]= temp_hess0;
      complex_double temp_hess1 = (hess1[0]+hess1[2])+(hess1[1]+hess1[3])*_Complex_I;
      hess[9*n+1]= temp_hess1;
      complex_double temp_hess2 = (hess2[0]+hess2[2])+(hess2[1]+hess2[3])*_Complex_I;
      hess[9*n+2]= temp_hess2;
      hess[9*n+3] = 0;
      complex_double temp_hess4 = (hess4[0]+hess4[2])+(hess4[1]+hess4[3])*_Complex_I;
      hess[9*n+4]= temp_hess4;
      complex_double temp_hess5 = (hess5[0]+hess5[2])+(hess5[1]+hess5[3])*_Complex_I;
      hess[9*n+5]= temp_hess5;
      hess[9*n+6] = 0;
      hess[9*n+7] = 0;
      complex_double temp_hess8 = (hess8[0]+hess8[2])+(hess8[1]+hess8[3])*_Complex_I;
      hess[9*n+8]= temp_hess8;
  }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;

  double dxInvdxInv = dxInv*dxInv;
  double dyInvdyInv = dyInv*dyInv;
  double dzInvdzInv = dzInv*dzInv;
  double dxInvdyInv = dxInv*dyInv;
  double dxInvdzInv = dxInv*dzInv;
  double dyInvdzInv = dyInv*dzInv;


 num_splines = spline->num_splines;

 for (int n=0; n<num_splines; n++) {
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

#else
void eval_multi_UBspline_3d_z(const multi_UBspline_3d_z *spline, double x, double y, double z, complex_double* restrict vals)
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
  
  double tpx[4];
  double tpy[4];
  double tpz[4];
  double a[4];
  double b[4];
  double c[4];
  double da[4];
  double db[4];
  double dc[4];
  double d2a[4];
  double d2b[4];
  double d2c[4];


  
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


  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;



  complex_double* restrict coefs = spline->coefs + ix*xs + iy*ys + iz*zs;
  int num_splines = spline->num_splines;
  for (int n=0; n<num_splines; n++, coefs++) {
      complex_double val = 0, grad0 = 0, grad1 = 0, grad2 = 0;
      complex_double hess0 = 0, hess1 = 0, hess2 = 0, hess4 = 0, hess5 = 0, hess8 = 0;

      for (int i=0; i<4; i++) {

	      double pre0 =   a[i] *   b[0];
	      double pre1 =  da[i] *   b[0];
	      double pre2 = d2a[i] *   b[0];
              double pre3 =   a[i] *  db[0];
	      complex_double coef0 = coefs[i*xs];
	      complex_double coef1 = coefs[i*xs + zs];
	      complex_double coef2 = coefs[i*xs + 2*zs];
	      complex_double coef3 = coefs[i*xs + 3*zs];
	      complex_double sum0 = c[0] * coef0 + c[1] * coef1 + c[2] * coef2 + c[3] * coef3;
	      complex_double sum1 = dc[0] * coef0 + dc[1] * coef1 + dc[2] * coef2 + dc[3] * coef3;

	      double pre01 =   a[i] *   b[1];
	      double pre11 =  da[i] *   b[1];
	      double pre21 = d2a[i] *   b[1];
              double pre31 =   a[i] *  db[1];
	      complex_double coef01 = coefs[i*xs + ys];
	      complex_double coef11 = coefs[i*xs + ys + zs];
	      complex_double coef21 = coefs[i*xs + ys + 2*zs];
	      complex_double coef31 = coefs[i*xs + ys + 3*zs];
	      complex_double sum01 = c[0] * coef01 + c[1] * coef11 + c[2] * coef21 + c[3] * coef31;
	      complex_double sum11 = dc[0] * coef01 + dc[1] * coef11 + dc[2] * coef21 + dc[3] * coef31;

	      double pre02 =   a[i] *   b[2];
	      double pre12 =  da[i] *   b[2];
	      double pre22 = d2a[i] *   b[2];
              double pre32 =   a[i] *  db[2];
	      complex_double coef02 = coefs[i*xs + 2*ys];
	      complex_double coef12 = coefs[i*xs + 2*ys + zs];
	      complex_double coef22 = coefs[i*xs + 2*ys + 2*zs];
	      complex_double coef32 = coefs[i*xs + 2*ys + 3*zs];
	      complex_double sum02 = c[0] * coef02 + c[1] * coef12 + c[2] * coef22 + c[3] * coef32;
	      complex_double sum12 = dc[0] * coef02 + dc[1] * coef12 + dc[2] * coef22 + dc[3] * coef32;

	      double pre03 =   a[i] *   b[3];
	      double pre13 =  da[i] *   b[3];
	      double pre23 = d2a[i] *   b[3];
              double pre33 =   a[i] *  db[3];
	      complex_double coef03 = coefs[i*xs + 3*ys];
	      complex_double coef13 = coefs[i*xs + 3*ys + zs];
	      complex_double coef23 = coefs[i*xs + 3*ys + 2*zs];
 	      complex_double coef33 = coefs[i*xs + 3*ys + 3*zs];
	      complex_double sum03 = c[0] * coef03 + c[1] * coef13 + c[2] * coef23 + c[3] * coef33;
	      complex_double sum13 = dc[0] * coef03 + dc[1] * coef13 + dc[2] * coef23 + dc[3] * coef33;

	      val   +=   pre0 * sum0 + pre01 * sum01 + pre02 * sum02 + pre03 * sum03;
	      hess8 += pre0 * (d2c[0] * coef0 +d2c[1] * coef1 + d2c[2] * coef2  + d2c[3] * coef3) + 
		  pre01 * (d2c[0] * coef01 +d2c[1] * coef11 + d2c[2] * coef21  + d2c[3] * coef31) +
		  pre02 * (d2c[0] * coef02 +d2c[1] * coef12 + d2c[2] * coef22  + d2c[3] * coef32) +
		  pre03 * (d2c[0] * coef03 +d2c[1] * coef13 + d2c[2] * coef23  + d2c[3] * coef33);
	      hess1 += (da[i] *  db[0]) * sum0 + (da[i] *  db[1]) * sum01 + (da[i] *  db[2]) * sum02 + (da[i] *  db[3]) * sum03;
	      hess4 += (a[i] * d2b[0]) * sum0 + (a[i] * d2b[1]) * sum01 + (a[i] * d2b[2]) * sum02 + (a[i] * d2b[3]) * sum03;
	      grad0 +=   pre1 * sum0 + pre11 * sum01 + pre12 * sum02 + pre13 * sum03;
	      grad1 +=   pre3 * sum0 + pre33 * sum03 + pre31 * sum01 + pre32 * sum02;
	      grad2 +=   pre0 * sum1 + pre01 * sum11 + pre02 * sum12 + pre03 * sum13;
	      hess2 +=   pre1 * sum1 + pre11 * sum11 + pre12 * sum12 + pre13 * sum13;
	      hess5 +=   pre3 * sum1 + pre31 * sum11 + pre32 * sum12 + pre33 * sum13;
	      hess0 +=   pre2 * sum0 + pre21 * sum01 + pre22 * sum02 + pre23 * sum03;

      }

      vals[n] = val;
      grads[3*n+0] = grad0;
      grads[3*n+1] = grad1;
      grads[3*n+2] = grad2;
      hess[9*n+0] = hess0;
      hess[9*n+1] = hess1;
      hess[9*n+2] = hess2;
      hess[9*n+3] = 0;
      hess[9*n+4] = hess4;
      hess[9*n+5] = hess5;
      hess[9*n+6] = 0;
      hess[9*n+7] = 0;
      hess[9*n+8] = hess8;

  }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv; 
 
 num_splines = spline->num_splines;
 
 for (int n=0; n<num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess[9*n+0] *= dxInv*dxInv;
    hess[9*n+4] *= dyInv*dyInv;
    hess[9*n+8] *= dzInv*dzInv;
    hess[9*n+1] *= dxInv*dyInv;
    hess[9*n+2] *= dxInv*dzInv;
    hess[9*n+5] *= dyInv*dzInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess[9*n+3] = hess[9*n+1];
    hess[9*n+6] = hess[9*n+2];
    hess[9*n+7] = hess[9*n+5];
  }
}

#endif
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
