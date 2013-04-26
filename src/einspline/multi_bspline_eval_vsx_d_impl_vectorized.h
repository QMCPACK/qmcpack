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

#ifndef MULTI_BSPLINE_EVAL_VSX_D_H
#define MULTI_BSPLINE_EVAL_VSX_D_H

#include <math.h>
#include <stdio.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"
#include <altivec.h>

extern const double* restrict   Ad;
extern const double* restrict  dAd;
extern const double* restrict d2Ad;

/************************************************************/
/* 1D double-precision, real evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_d (multi_UBspline_1d_d *spline,
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
eval_multi_UBspline_1d_d_vg (multi_UBspline_1d_d *spline,
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
eval_multi_UBspline_1d_d_vgl (multi_UBspline_1d_d *spline,
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
eval_multi_UBspline_1d_d_vgh (multi_UBspline_1d_d *spline,
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
eval_multi_UBspline_2d_d (multi_UBspline_2d_d *spline,
                          double x, double y,
                          double* restrict vals)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  double tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
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
    for (int j=0; j<4; j++)
    {
      double prefactor = a[i]*b[j];
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++)
        vals[n] += prefactor*coefs[n];
    }
}


void
eval_multi_UBspline_2d_d_vg (multi_UBspline_2d_d *spline,
                             double x, double y,
                             double* restrict vals,
                             double* restrict grads)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = grads[2*n+2] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      double ab = a[i]*b[j];
      double dab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++)
      {
        vals [n]     +=   ab   *coefs[n];
        grads[2*n+0] +=  dab[0]*coefs[n];
        grads[2*n+1] +=  dab[1]*coefs[n];
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
  }
}

void
eval_multi_UBspline_2d_d_vgl (multi_UBspline_2d_d *spline,
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    lapl2[2*n+0] = lapl2[2*n+1] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      double ab = a[i]*b[j];
      double dab[2], d2ab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i]*  b[j];
      d2ab[1] =   a[i]*d2b[j];
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++)
      {
        vals[n]      +=   ab   *coefs[n];
        grads[2*n+0] +=  dab[0]*coefs[n];
        grads[2*n+1] +=  dab[1]*coefs[n];
        lapl2[2*n+0] += d2ab[0]*coefs[n];
        lapl2[2*n+1] += d2ab[1]*coefs[n];
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    lapl2[2*n+0] *= dxInv*dxInv;
    lapl2[2*n+1] *= dyInv*dyInv;
    lapl[n] = lapl2[2*n+0] + lapl2[2*n+1];
  }
}




void
eval_multi_UBspline_2d_d_vgh (multi_UBspline_2d_d *spline,
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;
  tpx[1] = tx*tx;
  tpx[2] = tx;
  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;
  tpy[1] = ty*ty;
  tpy[2] = ty;
  tpy[3] = 1.0;
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    for (int i=0; i<4; i++)
      hess[4*n+i] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      double ab = a[i]*b[j];
      double dab[2], d2ab[3];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i] *   b[j];
      d2ab[1] =  da[i] *  db[j];
      d2ab[2] =   a[i] * d2b[j];
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++)
      {
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
  for (int n=0; n<spline->num_splines; n++)
  {
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
eval_multi_UBspline_3d_d (multi_UBspline_3d_d *spline,
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
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
  double* restrict coefs = spline->coefs;
  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  for (int n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;
  double *v = vals;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      vector double k0 = vec_splats(a[i]*b[j]*c[0]);
      vector double k1 = vec_splats(a[i]*b[j]*c[1]);
      vector double k2 = vec_splats(a[i]*b[j]*c[2]);
      vector double k3 = vec_splats(a[i]*b[j]*c[3]);
      vector double* restrict c0 = (vector double*)(spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+0)*zs));
      vector double* restrict c1 = (vector double*)(spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+1)*zs));
      vector double* restrict c2 = (vector double*)(spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+2)*zs));
      vector double* restrict c3 = (vector double*)(spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+3)*zs));
      for (int n=0; n<spline->num_splines/2; n+=2)
      {
        vector double v01 = *(vector double*)&(vals[2*n]);
        vector double v23 = *(vector double*)&(vals[2*n+1]);
        v01 = vec_madd (c0[n],   k0, v01);
        v23 = vec_madd (c0[n+1], k0, v23);
        v01 = vec_madd (c1[n],   k1, v01);
        v23 = vec_madd (c1[n+1], k1, v23);
        v01 = vec_madd (c2[n],   k2, v01);
        v23 = vec_madd (c2[n+1], k2, v23);
        v01 = vec_madd (c3[n],   k3, v01);
        v23 = vec_madd (c3[n+1], k3, v23);
        *(vector double*)&(vals[2*n])   = v01 ;
        *(vector double*)&(vals[2*n+2]) = v23 ;
        // vals[n+0] += k0*c0[n];
        // vals[n+1] += k0*c0[n+1];
        // vals[n+2] += k0*c0[n+2];
        // vals[n+3] += k0*c0[n+3];
        // vals[n+0] += k1*c1[n];
        // vals[n+1] += k1*c1[n+1];
        // vals[n+2] += k1*c1[n+2];
        // vals[n+3] += k1*c1[n+3];
        // vals[n+0] += k2*c2[n];
        // vals[n+1] += k2*c2[n+1];
        // vals[n+2] += k2*c2[n+2];
        // vals[n+3] += k2*c2[n+3];
        // vals[n+0] += k3*c3[n];
        // vals[n+1] += k3*c3[n+1];
        // vals[n+2] += k3*c3[n+2];
        // vals[n+3] += k3*c3[n+3];
        // vals[n+0] += k0*c0[n]   + k1*c1[n]   + k2*c2[n]   + k3*c3[n];
        // vals[n+1] += k0*c0[n+1] + k1*c1[n+1] + k2*c2[n+1] + k3*c3[n+1];
        // vals[n+2] += k0*c0[n+2] + k1*c1[n+2] + k2*c2[n+2] + k3*c3[n+2];
        // vals[n+3] += k0*c0[n+3] + k1*c1[n+3] + k2*c2[n+3] + k3*c3[n+3];
      }
    }
}


void
eval_multi_UBspline_3d_d_vg (multi_UBspline_3d_d *spline,
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4],
         da[4], db[4], dc[4];
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        double abc = a[i]*b[j]*c[k];
        double dabc[3];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];
        double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<spline->num_splines; n++)
        {
          vals[n]      +=   abc   *coefs[n];
          grads[3*n+0] +=  dabc[0]*coefs[n];
          grads[3*n+1] +=  dabc[1]*coefs[n];
          grads[3*n+2] +=  dabc[2]*coefs[n];
        }
      }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
  }
}



void
eval_multi_UBspline_3d_d_vgl (multi_UBspline_3d_d *spline,
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4],
         da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    lapl3[3*n+0] = lapl3[3*n+1] = lapl3[3*n+2] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        double abc = a[i]*b[j]*c[k];
        double dabc[3], d2abc[3];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];
        d2abc[0] = d2a[i]*  b[j]*  c[k];
        d2abc[1] =   a[i]*d2b[j]*  c[k];
        d2abc[2] =   a[i]*  b[j]*d2c[k];
        double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<spline->num_splines; n++)
        {
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
  for (int n=0; n<spline->num_splines; n++)
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
eval_multi_UBspline_3d_d_vgh (multi_UBspline_3d_d *spline,
                              double x, double y, double z,
                              double* restrict vals,
                              double* restrict grads,
                              double* restrict hess)
{
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4],
         da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }
  vector unsigned char i0 = (vector unsigned char)(0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
                            0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17);
  // vector unsigned char i1 = (vector unsigned char)(0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
  // 						   0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f);
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
// 	double abc = a[i]*b[j]*c[k];
// 	double dabc[3], d2abc[6];
// 	dabc[0] = da[i]* b[j]* c[k];
// 	dabc[1] =  a[i]*db[j]* c[k];
// 	dabc[2] =  a[i]* b[j]*dc[k];
// 	d2abc[0] = d2a[i]*  b[j]*  c[k];
// 	d2abc[1] =  da[i]* db[j]*  c[k];
// 	d2abc[2] =  da[i]*  b[j]* dc[k];
// 	d2abc[3] =   a[i]*d2b[j]*  c[k];
// 	d2abc[4] =   a[i]* db[j]* dc[k];
// 	d2abc[5] =   a[i]*  b[j]*d2c[k];
      vector double abc_v[40];
      abc_v[0] = vec_splats(  a[i] *  b[j] *  c[0]);
      abc_v[1] = vec_splats(  a[i] *  b[j] *  c[1]);
      abc_v[2] = vec_splats(  a[i] *  b[j] *  c[2]);
      abc_v[3] = vec_splats(  a[i] *  b[j] *  c[3]);
      abc_v[4] = vec_splats( da[i] *  b[j] *  c[0]);
      abc_v[5] = vec_splats( da[i] *  b[j] *  c[1]);
      abc_v[6] = vec_splats( da[i] *  b[j] *  c[2]);
      abc_v[7] = vec_splats( da[i] *  b[j] *  c[3]);
      abc_v[8] = vec_splats(  a[i] * db[j]  *  c[0]);
      abc_v[9] = vec_splats(  a[i] * db[j]  *  c[1]);
      abc_v[10] = vec_splats(  a[i] * db[j] *  c[2]);
      abc_v[11] = vec_splats(  a[i] * db[j] *  c[3]);
      abc_v[12] = vec_splats(  a[i] *  b[j] * dc[0]);
      abc_v[13] = vec_splats(  a[i] *  b[j] * dc[1]);
      abc_v[14] = vec_splats(  a[i] *  b[j] * dc[2]);
      abc_v[15] = vec_splats(  a[i] *  b[j] * dc[3]);
      abc_v[16] = vec_splats(d2a[i] *  b[j] *  c[0]);
      abc_v[17] = vec_splats(d2a[i] *  b[j] *  c[1]);
      abc_v[18] = vec_splats(d2a[i] *  b[j] *  c[2]);
      abc_v[19] = vec_splats(d2a[i] *  b[j] *  c[3]);
      abc_v[20] = vec_splats( da[i] * db[j] *  c[0]);
      abc_v[21] = vec_splats( da[i] * db[j] *  c[1]);
      abc_v[22] = vec_splats( da[i] * db[j] *  c[2]);
      abc_v[23] = vec_splats( da[i] * db[j] *  c[3]);
      abc_v[24] = vec_splats( da[i] *  b[j] * dc[0]);
      abc_v[25] = vec_splats( da[i] *  b[j] * dc[1]);
      abc_v[26] = vec_splats( da[i] *  b[j] * dc[2]);
      abc_v[27] = vec_splats( da[i] *  b[j] * dc[3]);
      abc_v[28] = vec_splats(  a[i] *d2b[j] *  c[0]);
      abc_v[29] = vec_splats(  a[i] *d2b[j] *  c[1]);
      abc_v[30] = vec_splats(  a[i] *d2b[j] *  c[2]);
      abc_v[31] = vec_splats(  a[i] *d2b[j] *  c[3]);
      abc_v[32] = vec_splats(  a[i] * db[j] * dc[0]);
      abc_v[33] = vec_splats(  a[i] * db[j] * dc[1]);
      abc_v[34] = vec_splats(  a[i] * db[j] * dc[2]);
      abc_v[35] = vec_splats(  a[i] * db[j] * dc[3]);
      abc_v[36] = vec_splats(  a[i] *  b[j] *d2c[0]);
      abc_v[37] = vec_splats(  a[i] *  b[j] *d2c[1]);
      abc_v[38] = vec_splats(  a[i] *  b[j] *d2c[2]);
      abc_v[39] = vec_splats(  a[i] *  b[j] *d2c[3]);
      double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz)*zs);
      for (int n=0; n<spline->num_splines; n+=2)
      {
        vector double c0 = *(vector double*)&(coefs[n     ]);//vec_xld2 (n     ,coefs);
        vector double c1 = *(vector double*)&(coefs[n+  zs]);//vec_xld2 (n+  zs,coefs);
        vector double c2 = *(vector double*)&(coefs[n+2*zs]);//vec_xld2 (n+2*zs,coefs);
        vector double c3 = *(vector double*)&(coefs[n+3*zs]);//vec_xld2 (n+3*zs,coefs);
        vector double v = *(vector double*)&(vals[n]);
        v           = vec_madd (c0, abc_v[0], v);
        v           = vec_madd (c1, abc_v[1], v);
        v           = vec_madd (c2, abc_v[2], v);
        v           = vec_madd (c3, abc_v[3], v);
        *(vector double*)&(vals[n]) = v;
        vector double ga = *(vector double*)&(grads[3*n+0]);
        vector double gb = *(vector double*)&(grads[3*n+2]);
        vector double gc = *(vector double*)&(grads[3*n+4]);
        // vector float a, b;
        // vector float c = vec_perm(a,b,i0);
        vector double g0 = vec_permi (ga, gb, 1); // 01
        vector double g1 = vec_permi (ga, gc, 2); // 10
        vector double g2 = vec_permi (gb, gc, 1); // 01
        g0 = vec_madd (c0, abc_v[ 4], g0);
        g0 = vec_madd (c1, abc_v[ 5], g0);
        g0 = vec_madd (c2, abc_v[ 6], g0);
        g0 = vec_madd (c3, abc_v[ 7], g0);
        g1 = vec_madd (c0, abc_v[ 8], g1);
        g1 = vec_madd (c1, abc_v[ 9], g1);
        g1 = vec_madd (c2, abc_v[10], g1);
        g1 = vec_madd (c3, abc_v[11], g1);
        g2 = vec_madd (c0, abc_v[12], g2);
        g2 = vec_madd (c1, abc_v[13], g2);
        g2 = vec_madd (c2, abc_v[14], g2);
        g2 = vec_madd (c3, abc_v[15], g2);
        *(vector double*)&(grads[3*n+0]) = vec_permi(g0,g1,0); // 00
        *(vector double*)&(grads[3*n+2]) = vec_permi(g2,g0,1); // 01
        *(vector double*)&(grads[3*n+4]) = vec_permi(g1,g2,3); // 11
        vector double h01 = *(vector double*)&(hess[9*n+0]);
        vector double h23 = *(vector double*)&(hess[9*n+2]);
        vector double h45 = *(vector double*)&(hess[9*n+4]);
        vector double h67 = *(vector double*)&(hess[9*n+6]);
        vector double h89 = *(vector double*)&(hess[9*n+8]);
        vector double h1011 = *(vector double*)&(hess[9*n+10]);
        vector double h1213 = *(vector double*)&(hess[9*n+12]);
        vector double h1415 = *(vector double*)&(hess[9*n+14]);
        vector double h1617 = *(vector double*)&(hess[9*n+16]);
        vector double h0 = vec_permi(h01,h89,  1);
        vector double h1 = vec_permi(h01,h1011,2);
        vector double h2 = vec_permi(h23,h1011,1);
        vector double h4 = vec_permi(h45,h1213,1);
        vector double h5 = vec_permi(h45,h1415,2);
        vector double h8 = vec_permi(h89,h1617,1);
        h0 = vec_madd (c0, abc_v[16], h0);
        h0 = vec_madd (c1, abc_v[17], h0);
        h0 = vec_madd (c2, abc_v[18], h0);
        h0 = vec_madd (c3, abc_v[19], h0);
        h1 = vec_madd (c0, abc_v[20], h1);
        h1 = vec_madd (c1, abc_v[21], h1);
        h1 = vec_madd (c2, abc_v[22], h1);
        h1 = vec_madd (c3, abc_v[23], h1);
        h2 = vec_madd (c0, abc_v[24], h2);
        h2 = vec_madd (c1, abc_v[25], h2);
        h2 = vec_madd (c2, abc_v[26], h2);
        h2 = vec_madd (c3, abc_v[27], h2);
        h4 = vec_madd (c0, abc_v[28], h4);
        h4 = vec_madd (c1, abc_v[29], h4);
        h4 = vec_madd (c2, abc_v[30], h4);
        h4 = vec_madd (c3, abc_v[31], h4);
        h5 = vec_madd (c0, abc_v[32], h5);
        h5 = vec_madd (c1, abc_v[33], h5);
        h5 = vec_madd (c2, abc_v[34], h5);
        h5 = vec_madd (c3, abc_v[35], h5);
        h8 = vec_madd (c0, abc_v[36], h8);
        h8 = vec_madd (c1, abc_v[37], h8);
        h8 = vec_madd (c2, abc_v[38], h8);
        h8 = vec_madd (c3, abc_v[39], h8);
        *(vector double*)&(hess[9*n+0])  = vec_permi(h0,h1,0);
        *(vector double*)&(hess[9*n+2])  = vec_permi(h2,h1,0);
        *(vector double*)&(hess[9*n+4])  = vec_permi(h4,h5,0);
        *(vector double*)&(hess[9*n+6])  = vec_permi(h2,h5,0);
        *(vector double*)&(hess[9*n+8])  = vec_permi(h8,h0,1);
        *(vector double*)&(hess[9*n+10]) = vec_permi(h1,h2,3);
        *(vector double*)&(hess[9*n+12]) = vec_permi(h1,h4,3);
        *(vector double*)&(hess[9*n+14]) = vec_permi(h5,h2,3);
        *(vector double*)&(hess[9*n+16]) = vec_permi(h5,h8,3);
        // vector double *g0 = (vector double*)&(grads[3*n+0]);
        // vector double *g1 = (vector double*)&(grads[3*n+1]);
        // vector double *g2 = (vector double*)&(grads[3*n+2]);
        // *g0 = vec_madd (c, *g0, abc_v[1]);
        // *g1 = vec_madd (c, *g1, abc_v[2]);
        // *g2 = vec_madd (c, *g2, abc_v[3]);
// 	  vector double* h0 = (vector double*)&(hess [9*n+0]);
// 	  vector double* h1 = (vector double*)&(hess [9*n+1]);
// 	  vector double* h2 = (vector double*)&(hess [9*n+2]);
// 	  vector double* h4 = (vector double*)&(hess [9*n+4]);
// 	  vector double* h5 = (vector double*)&(hess [9*n+5]);
// 	  vector double* h8 = (vector double*)&(hess [9*n+8]);
// 	  *h0  = vec_madd (c, *h0, abc_v[4]);
// 	  *h1  = vec_madd (c, *h1, abc_v[5]);
// 	  *h2  = vec_madd (c, *h2, abc_v[6]);
// 	  *h4  = vec_madd (c, *h4, abc_v[7]);
// 	  *h5  = vec_madd (c, *h5, abc_v[8]);
// 	  *h8  = vec_madd (c, *h8, abc_v[9]);
// 	  vals[n]      +=   abc   *coefs[n];
// 	  grads[3*n+0] +=  dabc[0]*coefs[n];
// 	  grads[3*n+1] +=  dabc[1]*coefs[n];
// 	  grads[3*n+2] +=  dabc[2]*coefs[n];
// 	  hess [9*n+0] += d2abc[0]*coefs[n];
// 	  hess [9*n+1] += d2abc[1]*coefs[n];
// 	  hess [9*n+2] += d2abc[2]*coefs[n];
// 	  hess [9*n+4] += d2abc[3]*coefs[n];
// 	  hess [9*n+5] += d2abc[4]*coefs[n];
// 	  hess [9*n+8] += d2abc[5]*coefs[n];
      }
    }
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  for (int n=0; n<spline->num_splines; n++)
  {
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
