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

#include <math.h>
#include <stdio.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"
#include "multi_bspline_eval_d.h"

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
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
  double tpx[4], a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  double* restrict coefs = spline->coefs;
  
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  double* restrict coefs0 = spline->coefs + ((ix  )*xs);
  double* restrict coefs1 = spline->coefs + ((ix+1)*xs);
  double* restrict coefs2 = spline->coefs + ((ix+2)*xs);
  double* restrict coefs3 = spline->coefs + ((ix+3)*xs);

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
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
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
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
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
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
  double dxInv = spline->x_grid.delta_inv;

  double* restrict coefs0 = spline->coefs + ((ix  )*xs);
  double* restrict coefs1 = spline->coefs + ((ix+1)*xs);
  double* restrict coefs2 = spline->coefs + ((ix+2)*xs);
  double* restrict coefs3 = spline->coefs + ((ix+3)*xs);

  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
    grads[n] = (da[0] * coefs0[n] + da[1] * coefs1[n] + da[2] * coefs2[n] + da[3] * coefs3[n])*dxInv;
    lapl[n]  = (d2a[0] * coefs0[n] + d2a[1] * coefs1[n] + d2a[2] * coefs2[n] + d2a[3] * coefs3[n])*dxInv*dxInv;
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


#ifdef BGQPX
#ifdef __xlC__
#include  <builtins.h>
#endif

void
eval_multi_UBspline_3d_d (const multi_UBspline_3d_d *spline,
			  double x, double y, double z,
			  double* restrict vals)
{
  // algorithm modified by Ye Luo, Mar. 12th 2015
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modf (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modf (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  
  double  a[4], b[4], c[4];

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

  vector4double vec_c0 = vec_splats(c[0]);
  vector4double vec_c1 = vec_splats(c[1]);
  vector4double vec_c2 = vec_splats(c[2]);
  vector4double vec_c3 = vec_splats(c[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      double pre00 =  a[i]*b[j];
      vector4double vec_pre00 = vec_splats(pre00);
      double* restrict coefs0 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      double* restrict coefs1 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs);
      double* restrict coefs2 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 2);
      double* restrict coefs3 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 3);

      int rest = spline->num_splines % 4;
      int n = 0;
      /* check alightment
      if ((long) &coefs0[n] & 0x1F ) printf("n0 is not 32­byte aligned\n");
      if ((long) &coefs1[n] & 0x1F ) printf("n1 is not 32­byte aligned\n");
      if ((long) &coefs2[n] & 0x1F ) printf("n2 is not 32­byte aligned\n");
      if ((long) &coefs3[n] & 0x1F ) printf("n3 is not 32­byte aligned\n";
      */
      for (int p=0; n<spline->num_splines-rest; n+=4, p+=32){
        vector4double vec_coef0, vec_coef1, vec_coef2, vec_coef3, vec_val;

        __dcbt(&coefs0[n+8]);
        __dcbt(&coefs1[n+8]);
        __dcbt(&coefs2[n+8]);
        __dcbt(&coefs3[n+8]);

        vec_coef0 = vec_ld(p, coefs0);
        vec_coef1 = vec_ld(p, coefs1);
        vec_coef2 = vec_ld(p, coefs2);
        vec_coef3 = vec_ld(p, coefs3);
        vec_val = vec_ld(p, vals);

        vec_coef0 = vec_mul(vec_c0, vec_coef0);
        vec_coef0 = vec_madd(vec_c1, vec_coef1, vec_coef0);
        vec_coef0 = vec_madd(vec_c2, vec_coef2, vec_coef0);
        vec_coef0 = vec_madd(vec_c3, vec_coef3, vec_coef0);
        vec_val = vec_madd(vec_pre00, vec_coef0, vec_val);
        vec_st(vec_val, p, vals);

      }
      for (; n<spline->num_splines; n++)
          vals[n] += pre00*(c[0]*coefs0[n] + c[1]*coefs1[n] + c[2]*coefs2[n] + c[3]*coefs3[n]);

    }
}

void
eval_multi_UBspline_3d_d_vgh (const multi_UBspline_3d_d *spline,
			      double x, double y, double z,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict hess) 
{
  // algorithm modified by Ye Luo, Mar. 12th 2015
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;

  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;

  //tx = modf(ux, &ipartx);  int ix = (int) ipartx;
  //ty = modf(uy, &iparty);  int iy = (int) iparty;
  //tz = modf(uz, &ipartz);  int iz = (int) ipartz;
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modf (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modf (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  
  double a[4], b[4], c[4],da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

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

  vector4double vec_c0 = vec_splats(c[0]);
  vector4double vec_c1 = vec_splats(c[1]);
  vector4double vec_c2 = vec_splats(c[2]);
  vector4double vec_c3 = vec_splats(c[3]);
  vector4double vec_dc0 = vec_splats(dc[0]);
  vector4double vec_dc1 = vec_splats(dc[1]);
  vector4double vec_dc2 = vec_splats(dc[2]);
  vector4double vec_dc3 = vec_splats(dc[3]);
  vector4double vec_d2c0 = vec_splats(d2c[0]);
  vector4double vec_d2c1 = vec_splats(d2c[1]);
  vector4double vec_d2c2 = vec_splats(d2c[2]);
  vector4double vec_d2c3 = vec_splats(d2c[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;

  }

  int rest = spline->num_splines % 4;
  int n = 0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      double* restrict coefs0 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      double* restrict coefs1 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs);
      double* restrict coefs2 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 2);
      double* restrict coefs3 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 3);

      double pre20 = d2a[i]*  b[j];
      double pre10 =  da[i]*  b[j];
      double pre00 =   a[i]*  b[j];
      double pre11 =  da[i]* db[j];
      double pre01 =   a[i]* db[j];
      double pre02 =   a[i]*d2b[j];

      vector4double vec_pre00 = vec_splats(pre00);
      vector4double vec_pre01 = vec_splats(pre01);
      vector4double vec_pre02 = vec_splats(pre02);
      vector4double vec_pre11 = vec_splats(pre11);
      vector4double vec_pre10 = vec_splats(pre10);
      vector4double vec_pre20 = vec_splats(pre20);

      n = 0;
      for (int hess_p=0, val_p=0; n<spline->num_splines-rest; n+=4, val_p+=32, hess_p+=288) {

        __dcbt(&coefs0[n+8]);
        __dcbt(&coefs1[n+8]);
        __dcbt(&coefs2[n+8]);
        __dcbt(&coefs3[n+8]);
        __dcbt(&vals  [n+8]);

        vector4double coef0 = vec_ld(0, &coefs0[n]);
        vector4double coef1 = vec_ld(0, &coefs1[n]);
        vector4double coef2 = vec_ld(0, &coefs2[n]);
        vector4double coef3 = vec_ld(0, &coefs3[n]);

        /* check alightment
        if ((long) &coefs1[n] & 0x1F ) fprintf(stderr,"n0 is not 32­byte aligned\n");
        if ((long) &coefs2[n] & 0x1F ) fprintf(stderr,"n1 is not 32­byte aligned\n");
        if ((long) &coefs3[n] & 0x1F ) fprintf(stderr,"n2 is not 32­byte aligned\n");
        if ((long) &coefs4[n] & 0x1F ) fprintf(stderr,"n3 is not 32­byte aligned\n");
        if ((long) hess & 0x1F )       fprintf(stderr,"hess is not 32­byte aligned\n");
        if ((long) vals & 0x1F )       fprintf(stderr,"vals is not 32­byte aligned\n");
        //if ((long) grads & 0x1F )    fprintf(stderr,"grads is not 32­byte aligned\n");
        */

        vector4double sum0, sum1, sum2;
        sum0 = vec_mul (vec_c0, coef0);
        sum0 = vec_madd(vec_c1, coef1, sum0);
        sum0 = vec_madd(vec_c2, coef2, sum0);
        sum0 = vec_madd(vec_c3, coef3, sum0);
        sum1 = vec_mul (vec_dc0, coef0);
        sum1 = vec_madd(vec_dc1, coef1, sum1);
        sum1 = vec_madd(vec_dc2, coef2, sum1);
        sum1 = vec_madd(vec_dc3, coef3, sum1);
        sum2 = vec_mul (vec_d2c0, coef0);
        sum2 = vec_madd(vec_d2c1, coef1, sum2);
        sum2 = vec_madd(vec_d2c2, coef2, sum2);
        sum2 = vec_madd(vec_d2c3, coef3, sum2);

        vector4double temp_vec;
        temp_vec = vec_ld(hess_p, hess);
        temp_vec = vec_madd(vec_pre20, sum0, temp_vec);
        vec_st(temp_vec, hess_p, hess);
        temp_vec = vec_ld(hess_p + 32, hess);
        temp_vec = vec_madd(vec_pre11, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 32, hess);
        temp_vec = vec_ld(hess_p + 64, hess);
        temp_vec = vec_madd(vec_pre10, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 64, hess);
        temp_vec = vec_ld(hess_p + 96, hess);
        temp_vec = vec_madd(vec_pre02, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 96, hess);
        temp_vec = vec_ld(hess_p + 128, hess);
        temp_vec = vec_madd(vec_pre01, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 128, hess);
        temp_vec = vec_ld(hess_p + 160, hess);
        temp_vec = vec_madd(vec_pre00, sum2, temp_vec);
        vec_st(temp_vec, hess_p + 160, hess);
        temp_vec = vec_ld(hess_p + 192, hess);
        temp_vec = vec_madd(vec_pre10, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 192, hess);
        temp_vec = vec_ld(hess_p + 224, hess);
        temp_vec = vec_madd(vec_pre01, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 224, hess);
        temp_vec = vec_ld(hess_p + 256, hess);
        temp_vec = vec_madd(vec_pre00, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 256, hess);
        temp_vec = vec_ld(val_p, vals);
        temp_vec = vec_madd(vec_pre00, sum0, temp_vec);
        vec_st(temp_vec, val_p, vals);
      }
      for (; n<spline->num_splines; n++) {
        double sum0 =   c[0] * coefs0[n] +   c[1] * coefs1[n] +   c[2] * coefs2[n] +   c[3] * coefs3[n];
        double sum1 =  dc[0] * coefs0[n] +  dc[1] * coefs1[n] +  dc[2] * coefs2[n] +  dc[3] * coefs3[n];
        double sum2 = d2c[0] * coefs0[n] + d2c[1] * coefs1[n] + d2c[2] * coefs2[n] + d2c[3] * coefs3[n];
        hess [9*n  ] += pre20 * sum0;
        hess [9*n+1] += pre11 * sum0;
        hess [9*n+2] += pre10 * sum1;
        hess [9*n+4] += pre02 * sum0;
        hess [9*n+5] += pre01 * sum1;
        hess [9*n+8] += pre00 * sum2;
        grads[3*n  ] += pre10 * sum0;
        grads[3*n+1] += pre01 * sum0;
        grads[3*n+2] += pre00 * sum1;
        vals[n]      += pre00 * sum0;
      }
    }

  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv; 

  n = 0;
  for(; n<spline->num_splines-rest; n+=4) {

    // reshuffle all contents in hess to hess and grads
    grads[3*n  ] = hess[9*n+24] * dxInv;
    grads[3*n+1] = hess[9*n+28] * dyInv;
    grads[3*n+2] = hess[9*n+32] * dzInv;
    grads[3*n+3] = hess[9*n+25] * dxInv;
    grads[3*n+4] = hess[9*n+29] * dyInv;
    grads[3*n+5] = hess[9*n+33] * dzInv;
    grads[3*n+6] = hess[9*n+26] * dxInv;
    grads[3*n+7] = hess[9*n+30] * dyInv;
    grads[3*n+8] = hess[9*n+34] * dzInv;
    grads[3*n+ 9] = hess[9*n+27] * dxInv;
    grads[3*n+10] = hess[9*n+31] * dyInv;
    grads[3*n+11] = hess[9*n+35] * dzInv;

                   hess[9*n+35] = hess[9*n+23] * dzInv*dzInv;
    hess[9*n+34] = hess[9*n+32] = hess[9*n+19] * dyInv*dzInv;
    hess[9*n+33] = hess[9*n+29] = hess[9*n+11] * dxInv*dzInv;
                   hess[9*n+31] = hess[9*n+15] * dyInv*dyInv;
    hess[9*n+30] = hess[9*n+28] = hess[9*n+ 7] * dxInv*dyInv;
                   hess[9*n+27] = hess[9*n+ 3] * dxInv*dxInv;
                   hess[9*n+26] = hess[9*n+22] * dzInv*dzInv;
    hess[9*n+25] = hess[9*n+23] = hess[9*n+18] * dyInv*dzInv;
                   hess[9*n+24] = hess[9*n+10] * dxInv*dzInv; // missing 20
                   hess[9*n+22] = hess[9*n+14] * dyInv*dyInv;
                   hess[9*n+18] = hess[9*n+ 2] * dxInv*dxInv;
                   hess[9*n+14] = hess[9*n+17] * dyInv*dzInv; // missing 16
                   hess[9*n+17] = hess[9*n+21] * dzInv*dzInv;
    hess[9*n+21] = hess[9*n+19] = hess[9*n+ 6] * dxInv*dyInv;
    hess[9*n+15] = hess[9*n+11] = hess[9*n+ 9] * dxInv*dzInv;
                   hess[9*n+10] = hess[9*n+ 5] * dxInv*dyInv; // mising 12
                   hess[9*n+ 9] = hess[9*n+ 1] * dxInv*dxInv;
    hess[9*n+ 3] = hess[9*n+ 1] = hess[9*n+ 4] * dxInv*dyInv;
    hess[9*n+ 2] = hess[9*n+ 6] = hess[9*n+ 8] * dxInv*dzInv;
                   hess[9*n+ 8] = hess[9*n+20] * dzInv*dzInv;
                   hess[9*n+20] = hess[9*n+24];
                   hess[9*n+ 4] = hess[9*n+12] * dyInv*dyInv;
    hess[9*n+ 7] = hess[9*n+ 5] = hess[9*n+16] * dyInv*dzInv; // missing 7
                   hess[9*n+12] = hess[9*n+10];
                   hess[9*n+16] = hess[9*n+14];
                   hess[9*n   ] = hess[9*n   ] * dxInv*dxInv;
                   hess[9*n+13] = hess[9*n+13] * dyInv*dyInv;
  }

  for (; n<spline->num_splines; n++) {

    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess [9*n+0] *= dxInv*dxInv;
    hess [9*n+1] *= dxInv*dyInv;
    hess [9*n+2] *= dxInv*dzInv;
    hess [9*n+4] *= dyInv*dyInv;
    hess [9*n+5] *= dyInv*dzInv;
    hess [9*n+8] *= dzInv*dzInv;
    // Copy hessian elements into lower half of 3x3 matrix
    hess [9*n+3] = hess[9*n+1];
    hess [9*n+6] = hess[9*n+2];
    hess [9*n+7] = hess[9*n+5];
  }
}

#else

void
eval_multi_UBspline_3d_d (const multi_UBspline_3d_d *spline,
			  double x, double y, double z,
			  double* restrict vals)
{
  // algorithm modified by Ye Luo, Mar. 12th 2015
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modf (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modf (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
   
  double  a[4], b[4], c[4];

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

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      double pre00 =  a[i]*b[j];
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      for(int n=0; n<spline->num_splines; n++)
          vals[n] += pre00*(c[0]*coefs[n] + c[1]*coefs[n+zs] + c[2]*coefs[n+2*zs] + c[3]*coefs[n+3*zs]);
    }
}

void
eval_multi_UBspline_3d_d_vgh (const multi_UBspline_3d_d *spline,
			      double x, double y, double z,
			      double* restrict vals,
			      double* restrict grads,
			      double* restrict hess) 
{
  // algorithm modified by Ye Luo, Mar. 12th 2015
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;

  //tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  //ty = modf (uy, &iparty);  int iy = (int) iparty;
  //tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modf (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modf (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);

  double a[4], b[4], c[4],da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  double* restrict coefs = spline->coefs;

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

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);

      double pre20 = d2a[i]*  b[j];
      double pre10 =  da[i]*  b[j];
      double pre00 =   a[i]*  b[j];
      double pre11 =  da[i]* db[j];
      double pre01 =   a[i]* db[j];
      double pre02 =   a[i]*d2b[j];

      for (int n=0; n<spline->num_splines; n++) {
        double sum0 =   c[0] * coefs[n] +   c[1] * coefs[n+zs] +   c[2] * coefs[n+2*zs] +   c[3] * coefs[n+3*zs];
        double sum1 =  dc[0] * coefs[n] +  dc[1] * coefs[n+zs] +  dc[2] * coefs[n+2*zs] +  dc[3] * coefs[n+3*zs];
        double sum2 = d2c[0] * coefs[n] + d2c[1] * coefs[n+zs] + d2c[2] * coefs[n+2*zs] + d2c[3] * coefs[n+3*zs];
        hess [9*n  ] += pre20 * sum0;
        hess [9*n+1] += pre11 * sum0;
        hess [9*n+2] += pre10 * sum1;
        hess [9*n+4] += pre02 * sum0;
        hess [9*n+5] += pre01 * sum1;
        hess [9*n+8] += pre00 * sum2;
        grads[3*n  ] += pre10 * sum0;
        grads[3*n+1] += pre01 * sum0;
        grads[3*n+2] += pre00 * sum1;
        vals[n]      += pre00 * sum0;
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
  }
}

#endif

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
  //tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  //ty = modf (uy, &iparty);  int iy = (int) iparty;
  //tz = modf (uz, &ipartz);  int iz = (int) ipartz;
  tx = modf (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modf (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modf (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);  

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
