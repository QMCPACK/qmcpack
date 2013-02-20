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

#ifndef MULTI_BSPLINE_EVAL_STD2_D_H
#define MULTI_BSPLINE_EVAL_STD2_D_H

#include <math.h>
#include <stdio.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"

extern const double* restrict   Ad;
extern const double* restrict  dAd;
extern const double* restrict d2Ad;
extern const double* restrict d3Ad;

/************************************************************/
/* 1D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_d (const multi_UBspline_1d_d *spline,
			  double x,
			  double* restrict vals)
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  double* restrict coefs = spline->coefs;
  
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  for (int n=0; n<spline->num_splines; n++) 
    vals[n]  = 0.0;

  for (int i=0; i<4; i++) {
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
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
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

  for (int n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }

  for (int i=0; i<4; i++) { 
    double* restrict coefs = spline->coefs + ((ix+i)*xs);
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
eval_multi_UBspline_1d_d_vgl (const multi_UBspline_1d_d *spline,
			      double x,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict lapl)	  
{
  x -= spline->x_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double ipartx, tx;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  
  double tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
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

  for (int n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }

  for (int i=0; i<4; i++) {      
    double* restrict coefs = spline->coefs + ((ix+i)*xs);
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
eval_multi_UBspline_1d_d_vgh (const multi_UBspline_1d_d *spline,
			      double x,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict hess)
{
  eval_multi_UBspline_1d_d_vgl (spline, x, vals, grads, hess);
}



/************************************************************/
/* 2D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_d (const multi_UBspline_2d_d *spline,
			  double x, double y,
			  double* restrict vals)
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
  double* restrict coefs = spline->coefs;

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
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) 
	vals[n] += prefactor*coefs[n];
    }
}


void
eval_multi_UBspline_2d_d_vg (const multi_UBspline_2d_d *spline,
			     double x, double y,
			     double* restrict vals,
			     double* restrict grads)
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
  double* restrict coefs = spline->coefs;

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
      
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
eval_multi_UBspline_2d_d_vgl (const multi_UBspline_2d_d *spline,
			      double x, double y,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict lapl)	  
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

  double lapl2[2*spline->num_splines];
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
      
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
eval_multi_UBspline_2d_d_vgh (const multi_UBspline_2d_d *spline,
			      double x, double y,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict hess)	  
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

      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
/* 3D double-precision, real evaulation functions           */
/************************************************************/
void
eval_multi_UBspline_3d_d (const multi_UBspline_3d_d *spline,
                          double x, double y, double z,
                          double* restrict vals)

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
  
  double a[4], b[4], c[4];
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


  double* restrict coefs = spline->coefs;


  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) 
      for (int k=0; k<4; k++) {
	double prefactor = a[i]*b[j]*c[k];
	double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	for (int n=0; n<spline->num_splines; n++) 
	  vals[n] += prefactor*coefs[n];
      }
}



void
eval_multi_UBspline_3d_d_vgh (const multi_UBspline_3d_d *spline,
			      double x, double y, double z,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict hess)	  
{

    double ux, uy, uz, prefactor, ipartx, iparty, ipartz, tx, ty, tz;
    double dxInv,dxInvdxInv,dyInv,dyInvdyInv,dzInv,dzInvdzInv,dxInvdyInv,dxInvdzInv,dyInvdzInv;
    double  a[4], b[4],c[4], da[4], db[4],dc[4], d2a[4], d2b[4],d2c[4]; 
    double dVal[64],dGrad0[64],dGrad1[64],dGrad2[64],dHess0[64],dHess1[64],dHess2[64],dHess4[64],dHess5[64],dHess8[64];
    double sVal,sGrad0,sGrad1,sGrad2,sHess0,sHess1,sHess2,sHess4,sHess5,sHess8;
    
    double *mod_coefs[64];
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

    M =  spline -> num_splines ;
    v = &vals[0];
    g = &grads[0];
    h = &hess[0];
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
       p = mod_coefs[i];

       for ( n = 0; n < M - rem; n = n + 4 ) 
{
          
           double  aVal0, aVal1, aVal2, aVal3;
           double  b0, b1, b2, b3;

           double  aGrad00, aGrad01, aGrad02;
           double  aGrad10, aGrad11, aGrad12;
           double  aGrad20, aGrad21, aGrad22;
           double  aGrad30, aGrad31, aGrad32;

           double  aHess00, aHess01, aHess02, aHess04,aHess05,aHess08;
           double  aHess10, aHess11, aHess12, aHess14,aHess15,aHess18;
           double  aHess20, aHess21, aHess22, aHess24,aHess25,aHess28;
           double  aHess30, aHess31, aHess32, aHess34,aHess35,aHess38;

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

          //Grad01
          aGrad00 = g[3*(n+0)+0];
          aGrad01 = g[3*(n+0)+1];
          aGrad02 = g[3*(n+0)+2];
          aGrad10 = g[3*(n+1)+0];

          aGrad00 = aGrad00 + sGrad0 * b0;
          aGrad01 = aGrad01 + sGrad1 * b0;
          aGrad02 = aGrad02 + sGrad2 * b0;
          aGrad10 = aGrad10 + sGrad0 * b1;

          g[3*(n+0)+0] = aGrad00;
          g[3*(n+0)+1] = aGrad01;
          g[3*(n+0)+2] = aGrad02;
          g[3*(n+1)+0] = aGrad10;

///////////Grad02

          aGrad11 = g[3*(n+1)+1];
          aGrad12 = g[3*(n+1)+2];
          aGrad20 = g[3*(n+2)+0];
          aGrad21 = g[3*(n+2)+1];

          aGrad11 = aGrad11 + sGrad1 * b1;
          aGrad12 = aGrad12 + sGrad2 * b1;
          aGrad20 = aGrad20 + sGrad0 * b2;
          aGrad21 = aGrad21 + sGrad1 * b2;

          g[3*(n+1)+1] = aGrad11;
          g[3*(n+1)+2] = aGrad12;
          g[3*(n+2)+0] = aGrad20;
          g[3*(n+2)+1] = aGrad01;



//////////////Grad03
          aGrad22 = g[3*(n+2)+2];
          aGrad30 = g[3*(n+3)+0];
          aGrad31 = g[3*(n+3)+1];
          aGrad32 = g[3*(n+3)+2];


          aGrad22 = aGrad22 + sGrad2 * b2;
          aGrad30 = aGrad30 + sGrad0 * b3;
          aGrad31 = aGrad31 + sGrad1 * b3;
          aGrad32 = aGrad32 + sGrad2 * b3;


          g[3*(n+2)+2] = aGrad22;
          g[3*(n+3)+0] = aGrad30;
          g[3*(n+3)+1] = aGrad31;
          g[3*(n+3)+2] = aGrad32;


//////////////////////////////////////






          //hess1
          aHess00 = h[9*(n+0)+0];
          aHess01 = h[9*(n+0)+1];
          aHess02 = h[9*(n+0)+2];
          aHess04 = h[9*(n+0)+4];

          aHess00 = aHess00 + sHess0 * b0;
          aHess01 = aHess01 + sHess1 * b0;
          aHess02 = aHess02 + sHess2 * b0;
          aHess04 = aHess04 + sHess4 * b0;
 
          h[9*(n+0)+0] = aHess00;
          h[9*(n+0)+1] = aHess01;
          h[9*(n+0)+2] = aHess02;
          h[9*(n+0)+4] = aHess04;

          
          //hess2
          aHess05 = h[9*(n+0)+5];
          aHess08 = h[9*(n+0)+8];
          aHess10 = h[9*(n+1)+0];
          aHess11 = h[9*(n+1)+1];

          aHess05 = aHess05 + sHess5 * b0;
          aHess08 = aHess08 + sHess8 * b0;
          aHess10 = aHess10 + sHess0 * b1;
          aHess11 = aHess11 + sHess1 * b1;
 
          h[9*(n+0)+5] = aHess05;
          h[9*(n+0)+8] = aHess08;
          h[9*(n+1)+0] = aHess10;
          h[9*(n+1)+1] = aHess11;

          //hess3
          aHess12 = h[9*(n+1)+2];
          aHess14 = h[9*(n+1)+4];
          aHess15 = h[9*(n+1)+5];
          aHess18 = h[9*(n+1)+8];

          aHess12 = aHess12 + sHess2 * b1;
          aHess14 = aHess14 + sHess4 * b1;
          aHess15 = aHess15 + sHess5 * b1;
          aHess18 = aHess18 + sHess8 * b1;
 
          h[9*(n+1)+2] = aHess12;
          h[9*(n+1)+4] = aHess14;
          h[9*(n+1)+5] = aHess15;
          h[9*(n+1)+8] = aHess18;

          //hess4
          aHess20 = h[9*(n+2)+0];
          aHess21 = h[9*(n+2)+1];
          aHess22 = h[9*(n+2)+2];
          aHess24 = h[9*(n+2)+4];

          aHess20 = aHess20 + sHess0 * b2;
          aHess21 = aHess21 + sHess1 * b2;
          aHess22 = aHess22 + sHess2 * b2;
          aHess24 = aHess24 + sHess4 * b2;
 
          h[9*(n+2)+0] = aHess20;
          h[9*(n+2)+1] = aHess21;
          h[9*(n+2)+2] = aHess22;
          h[9*(n+2)+4] = aHess24;

          
          //hess5
          aHess25 = h[9*(n+2)+5];
          aHess28 = h[9*(n+2)+8];
          aHess30 = h[9*(n+3)+0];
          aHess31 = h[9*(n+3)+1];

          aHess25 = aHess25 + sHess5 * b2;
          aHess28 = aHess28 + sHess8 * b2;
          aHess30 = aHess30 + sHess0 * b3;
          aHess31 = aHess31 + sHess1 * b3;
 
          h[9*(n+2)+5] = aHess25;
          h[9*(n+2)+8] = aHess28;
          h[9*(n+3)+0] = aHess30;
          h[9*(n+3)+1] = aHess31;

          //hess6
          aHess32 = h[9*(n+3)+2];
          aHess34 = h[9*(n+3)+4];
          aHess35 = h[9*(n+3)+5];
          aHess38 = h[9*(n+3)+8];

          aHess32 = aHess32 + sHess2 * b3;
          aHess34 = aHess34 + sHess4 * b3;
          aHess35 = aHess35 + sHess5 * b3;
          aHess38 = aHess38 + sHess8 * b3;
 
          h[9*(n+3)+2] = aHess32;
          h[9*(n+3)+4] = aHess34;
          h[9*(n+3)+5] = aHess35;
          h[9*(n+3)+8] = aHess38;

       }
       for ( k = n; k < M; k++ ){
           v[k] = v[k] + sVal * p[k];
           g[3*k+0] = g[3*k+0] + sGrad0 * p[k];
           g[3*k+1] = g[3*k+1] + sGrad1 * p[k];
           g[3*k+2] = g[3*k+2] + sGrad2 * p[k];
           h[9*k+0] = h[9*k+0] + sHess0 * p[k];
           h[9*k+1] = h[9*k+1] + sHess1 * p[k];
           h[9*k+2] = h[9*k+2] + sHess2 * p[k];
           h[9*k+4] = h[9*k+4] + sHess4 * p[k];
           h[9*k+5] = h[9*k+5] + sHess5 * p[k];
           h[9*k+8] = h[9*k+8] + sHess8 * p[k];
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
 
    for (int n=0; n<M; n++) {
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
eval_multi_UBspline_3d_d_vg (const multi_UBspline_3d_d *spline,
			     double x, double y, double z,
			     double* restrict vals,
			     double* restrict grads)
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
  double* restrict coefs = spline->coefs;

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

	double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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
eval_multi_UBspline_3d_d_vgl (const multi_UBspline_3d_d *spline,
			      double x, double y, double z,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict lapl)	  
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

  double lapl3[3*spline->num_splines];
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

	double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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
eval_multi_UBspline_3d_d_vghgh (const multi_UBspline_3d_d *spline,
               double x, double y, double z,
               double* restrict vals,
               double* restrict grads,
               double* restrict hess,
               double* restrict gradhess)    
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

   double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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
