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
#include "multi_bspline_eval_s.h"

using std::min;
using std::max;

extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;
extern const float* restrict d3Af;

/************************************************************/
/* 1D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_s (const multi_UBspline_1d_s *spline,
			  float x,
			  float* restrict vals)
{
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
  float tpx[4], a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;
  
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  float* restrict coefs0 = spline->coefs + ((ix  )*xs);
  float* restrict coefs1 = spline->coefs + ((ix+1)*xs);
  float* restrict coefs2 = spline->coefs + ((ix+2)*xs);
  float* restrict coefs3 = spline->coefs + ((ix+3)*xs);

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
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
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
  float tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
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

  for (int n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }

  for (int i=0; i<4; i++) { 
    float* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++) {
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
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  
  float tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
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
  float dxInv = spline->x_grid.delta_inv;

  float* restrict coefs0 = spline->coefs + ((ix  )*xs);
  float* restrict coefs1 = spline->coefs + ((ix+1)*xs);
  float* restrict coefs2 = spline->coefs + ((ix+2)*xs);
  float* restrict coefs3 = spline->coefs + ((ix+3)*xs);

  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = a[0] * coefs0[n] + a[1] * coefs1[n] + a[2] * coefs2[n] + a[3] * coefs3[n];
    grads[n] = (da[0] * coefs0[n] + da[1] * coefs1[n] + da[2] * coefs2[n] + da[3] * coefs3[n])*dxInv;
    lapl[n]  = (d2a[0] * coefs0[n] + d2a[1] * coefs1[n] + d2a[2] * coefs2[n] + d2a[3] * coefs3[n])*dxInv*dxInv;
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
/* 2D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_2d_s (const multi_UBspline_2d_s *spline,
			  double x, double y,
			  float* restrict vals)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      float prefactor = a[i]*b[j];
      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) 
	vals[n] += prefactor*coefs[n];
    }
}


void
eval_multi_UBspline_2d_s_vg (const multi_UBspline_2d_s *spline,
			     double x, double y,
			     float* restrict vals,
			     float* restrict grads)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = grads[2*n+2] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      float ab = a[i]*b[j];
      float dab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      
      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals [n]     +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
  }
}

void
eval_multi_UBspline_2d_s_vgl (const multi_UBspline_2d_s *spline,
			      double x, double y,
			      float* restrict vals,
			      float* restrict grads,
			      float* restrict lapl)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
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

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  float lapl2[2*spline->num_splines];
  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    lapl2[2*n+0] = lapl2[2*n+1] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) {
      float ab = a[i]*b[j];
      float dab[2], d2ab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i]*  b[j];
      d2ab[1] =   a[i]*d2b[j];
      
      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals[n]      +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
	lapl2[2*n+0] += d2ab[0]*coefs[n];
	lapl2[2*n+1] += d2ab[1]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    lapl2[2*n+0] *= dxInv*dxInv;
    lapl2[2*n+1] *= dyInv*dyInv;
    lapl[n] = lapl2[2*n+0] + lapl2[2*n+1];
  }
}




void
eval_multi_UBspline_2d_s_vgh (const multi_UBspline_2d_s *spline,
			      double x, double y,
			      float* restrict vals,
			      float* restrict grads,
			      float* restrict hess)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
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

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

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
      float ab = a[i]*b[j];
      float dab[2], d2ab[3];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i] *   b[j];
      d2ab[1] =  da[i] *  db[j];
      d2ab[2] =   a[i] * d2b[j];

      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++) {
	vals[n]      +=   ab   *coefs[n];
	grads[2*n+0] +=  dab[0]*coefs[n];
	grads[2*n+1] +=  dab[1]*coefs[n];
	hess [4*n+0] += d2ab[0]*coefs[n];
	hess [4*n+1] += d2ab[1]*coefs[n];
	hess [4*n+3] += d2ab[2]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
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
/* 3D double-precision, real evaulation functions        */
/************************************************************/

#ifdef BGQPX
#ifdef __xlC__
#include <builtins.h>
#endif

void
eval_multi_UBspline_3d_s (const multi_UBspline_3d_s *spline,
			  float x, float  y, float z,
			  float* restrict vals)
{
  // algorithm modified by Ye Luo, Mar. 12th 2015

  float ux, uy, uz, ipartx, iparty, ipartz, tx, ty, tz, a[4], b[4],c[4];
  intptr_t xs, ys, zs;

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;

  ux = x*spline->x_grid.delta_inv;
  uy = y*spline->y_grid.delta_inv;
  uz = z*spline->z_grid.delta_inv;
  
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);

  a[0] = ( ( Af[0]  * tx + Af[1] ) * tx + Af[2] ) * tx + Af[3];
  a[1] = ( ( Af[4]  * tx + Af[5] ) * tx + Af[6] ) * tx + Af[7];
  a[2] = ( ( Af[8]  * tx + Af[9] ) * tx + Af[10] ) * tx + Af[11];
  a[3] = ( ( Af[12] * tx + Af[13] ) * tx + Af[14] ) * tx + Af[15];

  b[0] = ( ( Af[0]  * ty + Af[1] ) * ty + Af[2] ) * ty + Af[3];
  b[1] = ( ( Af[4]  * ty + Af[5] ) * ty + Af[6] ) * ty + Af[7];
  b[2] = ( ( Af[8]  * ty + Af[9] ) * ty + Af[10] ) * ty + Af[11];
  b[3] = ( ( Af[12] * ty + Af[13] ) * ty + Af[14] ) * ty + Af[15];

  c[0] = ( ( Af[0]  * tz + Af[1] ) * tz + Af[2] ) * tz + Af[3];
  c[1] = ( ( Af[4]  * tz + Af[5] ) * tz + Af[6] ) * tz + Af[7];
  c[2] = ( ( Af[8]  * tz + Af[9] ) * tz + Af[10] ) * tz + Af[11];
  c[3] = ( ( Af[12] * tz + Af[13] ) * tz + Af[14] ) * tz + Af[15];

  vector4double vec_c0 = vec_splats((double)c[0]);
  vector4double vec_c1 = vec_splats((double)c[1]);
  vector4double vec_c2 = vec_splats((double)c[2]);
  vector4double vec_c3 = vec_splats((double)c[3]);

  xs = spline->x_stride;
  ys = spline->y_stride;
  zs = spline->z_stride;

  // std::cout << "x,y,z "<< ix << ", " << iy << ", " << iz << std::endl;

  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      float pre00 = a[i]*b[j];
      vector4double vec_pre00 = vec_splats((double)pre00);
      float* restrict coefs0 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      float* restrict coefs1 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs);
      float* restrict coefs2 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 2);
      float* restrict coefs3 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 3);
      int rest = spline->num_splines % 4;
      int n = 0;
      /* check alightment
      if ((long) &coefs0[n] & 0xF ) printf("n0 is not 16-byte aligned\n");
      if ((long) &coefs1[n] & 0xF ) printf("n1 is not 16-byte aligned\n");
      if ((long) &coefs2[n] & 0xF ) printf("n2 is not 16-byte aligned\n");
      if ((long) &coefs3[n] & 0xF ) printf("n3 is not 16-byte aligned\n";
      */
      for (int p=0; n<spline->num_splines-rest; n+=4, p+=16){
        vector4double vec_coef0, vec_coef1, vec_coef2, vec_coef3, vec_val;

        __dcbt(&coefs0[n+16]);
        __dcbt(&coefs1[n+16]);
        __dcbt(&coefs2[n+16]);
        __dcbt(&coefs3[n+16]);

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
eval_multi_UBspline_3d_s_vgh (const multi_UBspline_3d_s *spline,
			      float x, float y, float z,
			      float* restrict vals,
			      float* restrict grads,
			      float* restrict hess)
{
  // algorithm modified by Ye Luo, Mar. 12th 2015

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  //tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  //ty = modff (uy, &iparty);  int iy = (int) iparty;
  //tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  
  float a[4], b[4], c[4],da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];

  a[0] = ( ( Af[0]  * tx + Af[1] ) * tx + Af[2] ) * tx + Af[3]; 
  a[1] = ( ( Af[4]  * tx + Af[5] ) * tx + Af[6] ) * tx + Af[7]; 
  a[2] = ( ( Af[8]  * tx + Af[9] ) * tx + Af[10] ) * tx + Af[11]; 
  a[3] = ( ( Af[12] * tx + Af[13] ) * tx + Af[14] ) * tx + Af[15]; 
  da[0] = ( ( dAf[0]  * tx + dAf[1] ) * tx + dAf[2] ) * tx + dAf[3]; 
  da[1] = ( ( dAf[4]  * tx + dAf[5] ) * tx + dAf[6] ) * tx + dAf[7]; 
  da[2] = ( ( dAf[8]  * tx + dAf[9] ) * tx + dAf[10] ) * tx + dAf[11]; 
  da[3] = ( ( dAf[12] * tx + dAf[13] ) * tx + dAf[14] ) * tx + dAf[15]; 
  d2a[0] = ( ( d2Af[0]  * tx + d2Af[1] ) * tx + d2Af[2] ) * tx + d2Af[3]; 
  d2a[1] = ( ( d2Af[4]  * tx + d2Af[5] ) * tx + d2Af[6] ) * tx + d2Af[7]; 
  d2a[2] = ( ( d2Af[8]  * tx + d2Af[9] ) * tx + d2Af[10] ) * tx + d2Af[11]; 
  d2a[3] = ( ( d2Af[12] * tx + d2Af[13] ) * tx + d2Af[14] ) * tx + d2Af[15]; 

  
  b[0] = ( ( Af[0]  * ty + Af[1] ) * ty + Af[2] ) * ty + Af[3]; 
  b[1] = ( ( Af[4]  * ty + Af[5] ) * ty + Af[6] ) * ty + Af[7]; 
  b[2] = ( ( Af[8]  * ty + Af[9] ) * ty + Af[10] ) * ty + Af[11]; 
  b[3] = ( ( Af[12] * ty + Af[13] ) * ty + Af[14] ) * ty + Af[15]; 
  db[0] = ( ( dAf[0]  * ty + dAf[1] ) * ty + dAf[2] ) * ty + dAf[3]; 
  db[1] = ( ( dAf[4]  * ty + dAf[5] ) * ty + dAf[6] ) * ty + dAf[7]; 
  db[2] = ( ( dAf[8]  * ty + dAf[9] ) * ty + dAf[10] ) * ty + dAf[11]; 
  db[3] = ( ( dAf[12] * ty + dAf[13] ) * ty + dAf[14] ) * ty + dAf[15]; 
  d2b[0] = ( ( d2Af[0]  * ty + d2Af[1] ) * ty + d2Af[2] ) * ty + d2Af[3]; 
  d2b[1] = ( ( d2Af[4]  * ty + d2Af[5] ) * ty + d2Af[6] ) * ty + d2Af[7]; 
  d2b[2] = ( ( d2Af[8]  * ty + d2Af[9] ) * ty + d2Af[10] ) * ty + d2Af[11]; 
  d2b[3] = ( ( d2Af[12] * ty + d2Af[13] ) * ty + d2Af[14] ) * ty + d2Af[15]; 
 
 

  c[0] = ( ( Af[0]  * tz + Af[1] ) * tz + Af[2] ) * tz + Af[3]; 
  c[1] = ( ( Af[4]  * tz + Af[5] ) * tz + Af[6] ) * tz + Af[7]; 
  c[2] = ( ( Af[8]  * tz + Af[9] ) * tz + Af[10] ) * tz + Af[11]; 
  c[3] = ( ( Af[12] * tz + Af[13] ) * tz + Af[14] ) * tz + Af[15]; 
  dc[0] = ( ( dAf[0]  * tz + dAf[1] ) * tz + dAf[2] ) * tz + dAf[3]; 
  dc[1] = ( ( dAf[4]  * tz + dAf[5] ) * tz + dAf[6] ) * tz + dAf[7]; 
  dc[2] = ( ( dAf[8]  * tz + dAf[9] ) * tz + dAf[10] ) * tz + dAf[11]; 
  dc[3] = ( ( dAf[12] * tz + dAf[13] ) * tz + dAf[14] ) * tz + dAf[15]; 
  d2c[0] = ( ( d2Af[0]  * tz + d2Af[1] ) * tz + d2Af[2] ) * tz + d2Af[3]; 
  d2c[1] = ( ( d2Af[4]  * tz + d2Af[5] ) * tz + d2Af[6] ) * tz + d2Af[7]; 
  d2c[2] = ( ( d2Af[8]  * tz + d2Af[9] ) * tz + d2Af[10] ) * tz + d2Af[11]; 
  d2c[3] = ( ( d2Af[12] * tz + d2Af[13] ) * tz + d2Af[14] ) * tz + d2Af[15]; 

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
      float* restrict coefs0 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      float* restrict coefs1 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs);
      float* restrict coefs2 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 2);
      float* restrict coefs3 = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs + zs * 3);

      float pre20 = d2a[i]*  b[j];
      float pre10 =  da[i]*  b[j];
      float pre00 =   a[i]*  b[j];
      float pre11 =  da[i]* db[j];
      float pre01 =   a[i]* db[j];
      float pre02 =   a[i]*d2b[j];

      vector4double vec_pre00 = vec_splats(pre00);
      vector4double vec_pre01 = vec_splats(pre01);
      vector4double vec_pre02 = vec_splats(pre02);
      vector4double vec_pre11 = vec_splats(pre11);
      vector4double vec_pre10 = vec_splats(pre10);
      vector4double vec_pre20 = vec_splats(pre20);

      n = 0;
      for (int hess_p=0, val_p=0; n<spline->num_splines-rest; n+=4, val_p+=16, hess_p+=144) {

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
        if ((long) &coefs0[n] & 0xF ) fprintf(stderr, "n0 is not 16-byte aligned\n");
        if ((long) &coefs1[n] & 0xF ) fprintf(stderr, "n1 is not 16-byte aligned\n");
        if ((long) &coefs2[n] & 0xF ) fprintf(stderr, "n2 is not 16-byte aligned\n");
        if ((long) &coefs3[n] & 0xF ) fprintf(stderr, "n3 is not 16-byte aligned\n");
        if ((long) hess & 0xF ) fprintf(stderr, "hess is not 16-byte aligned\n");
        if ((long) vals & 0xF ) fprintf(stderr, "vals is not 16-byte aligned\n");
        // if ((long) grads & 0xF ) fprintf(stderr,"grads is not 16-byte aligned\n");
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
        temp_vec = vec_ld(hess_p + 16, hess);
        temp_vec = vec_madd(vec_pre11, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 16, hess);
        temp_vec = vec_ld(hess_p + 32, hess);
        temp_vec = vec_madd(vec_pre10, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 32, hess);
        temp_vec = vec_ld(hess_p + 48, hess);
        temp_vec = vec_madd(vec_pre02, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 48, hess);
        temp_vec = vec_ld(hess_p + 64, hess);
        temp_vec = vec_madd(vec_pre01, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 64, hess);
        temp_vec = vec_ld(hess_p + 80, hess);
        temp_vec = vec_madd(vec_pre00, sum2, temp_vec);
        vec_st(temp_vec, hess_p + 80, hess);
        temp_vec = vec_ld(hess_p + 96, hess);
        temp_vec = vec_madd(vec_pre10, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 96, hess);
        temp_vec = vec_ld(hess_p + 112, hess);
        temp_vec = vec_madd(vec_pre01, sum0, temp_vec);
        vec_st(temp_vec, hess_p + 112, hess);
        temp_vec = vec_ld(hess_p + 128, hess);
        temp_vec = vec_madd(vec_pre00, sum1, temp_vec);
        vec_st(temp_vec, hess_p + 128, hess);
        temp_vec = vec_ld(val_p, vals);
        temp_vec = vec_madd(vec_pre00, sum0, temp_vec);
        vec_st(temp_vec, val_p, vals);
      }
      for (; n<spline->num_splines; n++) {
        float sum0 =   c[0] * coefs0[n] +   c[1] * coefs1[n] +   c[2] * coefs2[n] +   c[3] * coefs3[n];
        float sum1 =  dc[0] * coefs0[n] +  dc[1] * coefs1[n] +  dc[2] * coefs2[n] +  dc[3] * coefs3[n];
        float sum2 = d2c[0] * coefs0[n] + d2c[1] * coefs1[n] + d2c[2] * coefs2[n] + d2c[3] * coefs3[n];
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

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 

  n = 0;
  for (; n<spline->num_splines-rest; n+=4) {

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
eval_multi_UBspline_3d_s (const multi_UBspline_3d_s *spline,
			  float x, float  y, float z,
			  float* restrict vals)
{
  // algorithm modified by Ye Luo, Mar. 12th 2015

 float ux, uy, uz, ipartx, iparty, ipartz, tx, ty, tz, a[4], b[4],c[4];
 intptr_t xs, ys, zs;


  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;

  ux = x*spline->x_grid.delta_inv;
  uy = y*spline->y_grid.delta_inv;
  uz = z*spline->z_grid.delta_inv;
  
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);

  a[0] = ( ( Af[0]  * tx + Af[1] ) * tx + Af[2] ) * tx + Af[3];
  a[1] = ( ( Af[4]  * tx + Af[5] ) * tx + Af[6] ) * tx + Af[7];
  a[2] = ( ( Af[8]  * tx + Af[9] ) * tx + Af[10] ) * tx + Af[11];
  a[3] = ( ( Af[12] * tx + Af[13] ) * tx + Af[14] ) * tx + Af[15];

  b[0] = ( ( Af[0]  * ty + Af[1] ) * ty + Af[2] ) * ty + Af[3];
  b[1] = ( ( Af[4]  * ty + Af[5] ) * ty + Af[6] ) * ty + Af[7];
  b[2] = ( ( Af[8]  * ty + Af[9] ) * ty + Af[10] ) * ty + Af[11];
  b[3] = ( ( Af[12] * ty + Af[13] ) * ty + Af[14] ) * ty + Af[15];

  c[0] = ( ( Af[0]  * tz + Af[1] ) * tz + Af[2] ) * tz + Af[3];
  c[1] = ( ( Af[4]  * tz + Af[5] ) * tz + Af[6] ) * tz + Af[7];
  c[2] = ( ( Af[8]  * tz + Af[9] ) * tz + Af[10] ) * tz + Af[11];
  c[3] = ( ( Af[12] * tz + Af[13] ) * tz + Af[14] ) * tz + Af[15];

  xs = spline->x_stride;
  ys = spline->y_stride;
  zs = spline->z_stride;


  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++){
      float pre00 =  a[i]*b[j];
      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);
      for(int n=0; n<spline->num_splines; n++)
	  vals[n] += pre00*(c[0]*coefs[n] + c[1]*coefs[n+zs] + c[2]*coefs[n+2*zs] + c[3]*coefs[n+3*zs]);
    }

}

void
eval_multi_UBspline_3d_s_vgh (const multi_UBspline_3d_s *spline,
			      float x, float y, float z,
			      float* restrict vals,
			      float* restrict grads,
			      float* restrict hess) 
{
  // algorithm modified by Ye Luo, Mar. 12th 2015

  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  //tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  //ty = modff (uy, &iparty);  int iy = (int) iparty;
  //tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  
  float a[4], b[4], c[4],da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  float* restrict coefs = spline->coefs;
    a[0] = ( ( Af[0]  * tx + Af[1] ) * tx + Af[2] ) * tx + Af[3]; 
    a[1] = ( ( Af[4]  * tx + Af[5] ) * tx + Af[6] ) * tx + Af[7]; 
    a[2] = ( ( Af[8]  * tx + Af[9] ) * tx + Af[10] ) * tx + Af[11]; 
    a[3] = ( ( Af[12] * tx + Af[13] ) * tx + Af[14] ) * tx + Af[15]; 
    da[0] = ( ( dAf[0]  * tx + dAf[1] ) * tx + dAf[2] ) * tx + dAf[3]; 
    da[1] = ( ( dAf[4]  * tx + dAf[5] ) * tx + dAf[6] ) * tx + dAf[7]; 
    da[2] = ( ( dAf[8]  * tx + dAf[9] ) * tx + dAf[10] ) * tx + dAf[11]; 
    da[3] = ( ( dAf[12] * tx + dAf[13] ) * tx + dAf[14] ) * tx + dAf[15]; 
    d2a[0] = ( ( d2Af[0]  * tx + d2Af[1] ) * tx + d2Af[2] ) * tx + d2Af[3]; 
    d2a[1] = ( ( d2Af[4]  * tx + d2Af[5] ) * tx + d2Af[6] ) * tx + d2Af[7]; 
    d2a[2] = ( ( d2Af[8]  * tx + d2Af[9] ) * tx + d2Af[10] ) * tx + d2Af[11]; 
    d2a[3] = ( ( d2Af[12] * tx + d2Af[13] ) * tx + d2Af[14] ) * tx + d2Af[15]; 

    
    b[0] = ( ( Af[0]  * ty + Af[1] ) * ty + Af[2] ) * ty + Af[3]; 
    b[1] = ( ( Af[4]  * ty + Af[5] ) * ty + Af[6] ) * ty + Af[7]; 
    b[2] = ( ( Af[8]  * ty + Af[9] ) * ty + Af[10] ) * ty + Af[11]; 
    b[3] = ( ( Af[12] * ty + Af[13] ) * ty + Af[14] ) * ty + Af[15]; 
    db[0] = ( ( dAf[0]  * ty + dAf[1] ) * ty + dAf[2] ) * ty + dAf[3]; 
    db[1] = ( ( dAf[4]  * ty + dAf[5] ) * ty + dAf[6] ) * ty + dAf[7]; 
    db[2] = ( ( dAf[8]  * ty + dAf[9] ) * ty + dAf[10] ) * ty + dAf[11]; 
    db[3] = ( ( dAf[12] * ty + dAf[13] ) * ty + dAf[14] ) * ty + dAf[15]; 
    d2b[0] = ( ( d2Af[0]  * ty + d2Af[1] ) * ty + d2Af[2] ) * ty + d2Af[3]; 
    d2b[1] = ( ( d2Af[4]  * ty + d2Af[5] ) * ty + d2Af[6] ) * ty + d2Af[7]; 
    d2b[2] = ( ( d2Af[8]  * ty + d2Af[9] ) * ty + d2Af[10] ) * ty + d2Af[11]; 
    d2b[3] = ( ( d2Af[12] * ty + d2Af[13] ) * ty + d2Af[14] ) * ty + d2Af[15]; 
   
   
 
    c[0] = ( ( Af[0]  * tz + Af[1] ) * tz + Af[2] ) * tz + Af[3]; 
    c[1] = ( ( Af[4]  * tz + Af[5] ) * tz + Af[6] ) * tz + Af[7]; 
    c[2] = ( ( Af[8]  * tz + Af[9] ) * tz + Af[10] ) * tz + Af[11]; 
    c[3] = ( ( Af[12] * tz + Af[13] ) * tz + Af[14] ) * tz + Af[15]; 
    dc[0] = ( ( dAf[0]  * tz + dAf[1] ) * tz + dAf[2] ) * tz + dAf[3]; 
    dc[1] = ( ( dAf[4]  * tz + dAf[5] ) * tz + dAf[6] ) * tz + dAf[7]; 
    dc[2] = ( ( dAf[8]  * tz + dAf[9] ) * tz + dAf[10] ) * tz + dAf[11]; 
    dc[3] = ( ( dAf[12] * tz + dAf[13] ) * tz + dAf[14] ) * tz + dAf[15]; 
    d2c[0] = ( ( d2Af[0]  * tz + d2Af[1] ) * tz + d2Af[2] ) * tz + d2Af[3]; 
    d2c[1] = ( ( d2Af[4]  * tz + d2Af[5] ) * tz + d2Af[6] ) * tz + d2Af[7]; 
    d2c[2] = ( ( d2Af[8]  * tz + d2Af[9] ) * tz + d2Af[10] ) * tz + d2Af[11]; 
    d2c[3] = ( ( d2Af[12] * tz + d2Af[13] ) * tz + d2Af[14] ) * tz + d2Af[15]; 


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
      float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + iz*zs);

      float pre20 = d2a[i]*  b[j];
      float pre10 =  da[i]*  b[j];
      float pre00 =   a[i]*  b[j];
      float pre11 =  da[i]* db[j];
      float pre01 =   a[i]* db[j];
      float pre02 =   a[i]*d2b[j];

      for (int n=0; n<spline->num_splines; n++) {
        float sum0 =   c[0] * coefs[n] +   c[1] * coefs[n+zs] +   c[2] * coefs[n+2*zs] +   c[3] * coefs[n+3*zs];
        float sum1 =  dc[0] * coefs[n] +  dc[1] * coefs[n+zs] +  dc[2] * coefs[n+2*zs] +  dc[3] * coefs[n+3*zs];
        float sum2 = d2c[0] * coefs[n] + d2c[1] * coefs[n+zs] + d2c[2] * coefs[n+2*zs] + d2c[3] * coefs[n+3*zs];
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

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 

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
eval_multi_UBspline_3d_s_vg (const multi_UBspline_3d_s *spline,
			     float x, float y, float z,// double x, double y, double z,
			     float* restrict vals,
			     float* restrict grads)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);

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
	float abc = a[i]*b[j]*c[k];
	float dabc[3];
	dabc[0] = da[i]* b[j]* c[k];
	dabc[1] =  a[i]*db[j]* c[k];
	dabc[2] =  a[i]* b[j]*dc[k];

	float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
	for (int n=0; n<spline->num_splines; n++) {
	  vals[n]      +=   abc   *coefs[n];
	  grads[3*n+0] +=  dabc[0]*coefs[n];
	  grads[3*n+1] +=  dabc[1]*coefs[n];
	  grads[3*n+2] +=  dabc[2]*coefs[n];
	}
      }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
  for (int n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
  }
}



void
eval_multi_UBspline_3d_s_vgl (const multi_UBspline_3d_s *spline,
			      float x, float y, float z,// double x, double y, double z,
			      float* restrict vals,
			      float* restrict grads,
			      float* restrict lapl)	  
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  //tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  //ty = modff (uy, &iparty);  int iy = (int) iparty;
  //tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  tx = modff (ux, &ipartx);  int ix = std::min(std::max(0,(int) ipartx),spline->x_grid.num-1);
  ty = modff (uy, &iparty);  int iy = std::min(std::max(0,(int) iparty),spline->y_grid.num-1);
  tz = modff (uz, &ipartz);  int iz = std::min(std::max(0,(int) ipartz),spline->z_grid.num-1);
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
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

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 0]*tpz[0] + d2Af[ 1]*tpz[1] + d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 4]*tpz[0] + d2Af[ 5]*tpz[1] + d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[ 8]*tpz[0] + d2Af[ 9]*tpz[1] + d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[12]*tpz[0] + d2Af[13]*tpz[1] + d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  float lapl3[3*spline->num_splines];
  for (int n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    lapl3[3*n+0] = lapl3[3*n+1] = lapl3[3*n+2] = 0.0;
  }

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++) 
      for (int k=0; k<4; k++) {
	float abc = a[i]*b[j]*c[k];
	float dabc[3], d2abc[3];
	dabc[0] = da[i]* b[j]* c[k];
	dabc[1] =  a[i]*db[j]* c[k];
	dabc[2] =  a[i]* b[j]*dc[k];
	d2abc[0] = d2a[i]*  b[j]*  c[k];
	d2abc[1] =   a[i]*d2b[j]*  c[k];
	d2abc[2] =   a[i]*  b[j]*d2c[k];

	float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
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
eval_multi_UBspline_3d_s_vghgh (const multi_UBspline_3d_s *spline,
               float x, float y, float z,// double x, double y, double z,
               float* restrict vals,
               float* restrict grads,
               float* restrict hess,
               float* restrict gradhess)    
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4],
    d3a[4], d3b[4], d3c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
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
  d3a[0] = (/*d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] +*/ d3Af[ 3]*tpx[3]);
  d3a[1] = (/*d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] +*/ d3Af[ 7]*tpx[3]);
  d3a[2] = (/*d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] +*/ d3Af[11]*tpx[3]);
  d3a[3] = (/*d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] +*/ d3Af[15]*tpx[3]);
  
  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  d3b[0] = (/*d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] +*/ d3Af[ 3]*tpy[3]);
  d3b[1] = (/*d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] +*/ d3Af[ 7]*tpy[3]);
  d3b[2] = (/*d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] +*/ d3Af[11]*tpy[3]);
  d3b[3] = (/*d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] +*/ d3Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 0]*tpz[0] + d2Af[ 1]*tpz[1] + d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 4]*tpz[0] + d2Af[ 5]*tpz[1] + d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[ 8]*tpz[0] + d2Af[ 9]*tpz[1] + d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[12]*tpz[0] + d2Af[13]*tpz[1] + d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);
  d3c[0] = (/*d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] +*/ d3Af[ 3]*tpz[3]);
  d3c[1] = (/*d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] +*/ d3Af[ 7]*tpz[3]);
  d3c[2] = (/*d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] +*/ d3Af[11]*tpz[3]);
  d3c[3] = (/*d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] +*/ d3Af[15]*tpz[3]);

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
   float abc = a[i]*b[j]*c[k];
   float dabc[3], d2abc[6], d3abc[10];
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

   float* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
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
