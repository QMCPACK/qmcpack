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

#ifndef MULTI_BSPLINE_EVAL_STD_Z_IMPL_H
#define MULTI_BSPLINE_EVAL_STD_Z_IMPL_H

#include <math.h>
#include <stdio.h>
#include "bspline_base.h"
#include "multi_bspline_structs.h"

extern const double* restrict   Ad;
extern const double* restrict  dAd;
extern const double* restrict d2Ad;

/************************************************************/
/* 1D double-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_1d_z (multi_UBspline_1d_z *spline,
                          double x,
                          complex_double* restrict vals)
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
  complex_double* restrict coefs = spline->coefs;
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  intptr_t xs = spline->x_stride;
  for (int n=0; n<spline->num_splines; n++)
    vals[n]  = 0.0;
  for (int i=0; i<4; i++)
  {
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
    for (int n=0; n<spline->num_splines; n++)
      vals[n]  +=   a[i] * coefs[n];
  }
}



void
eval_multi_UBspline_1d_z_vg (multi_UBspline_1d_z *spline,
                             double x,
                             complex_double* restrict vals,
                             complex_double* restrict grads)
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
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
eval_multi_UBspline_1d_z_vgl (multi_UBspline_1d_z *spline,
                              double x,
                              complex_double* restrict vals,
                              complex_double* restrict grads,
                              complex_double* restrict lapl)
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }
  for (int i=0; i<4; i++)
  {
    complex_double* restrict coefs = spline->coefs + ((ix+i)*xs);
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
eval_multi_UBspline_1d_z_vgh (multi_UBspline_1d_z *spline,
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
eval_multi_UBspline_2d_z (multi_UBspline_2d_z *spline,
                          double x, double y,
                          complex_double* restrict vals)
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
    for (int j=0; j<4; j++)
    {
      double prefactor = a[i]*b[j];
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (int n=0; n<spline->num_splines; n++)
        vals[n] += prefactor*coefs[n];
    }
}


void
eval_multi_UBspline_2d_z_vg (multi_UBspline_2d_z *spline,
                             double x, double y,
                             complex_double* restrict vals,
                             complex_double* restrict grads)
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
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
eval_multi_UBspline_2d_z_vgl (multi_UBspline_2d_z *spline,
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
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
eval_multi_UBspline_2d_z_vgh (multi_UBspline_2d_z *spline,
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
      complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
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
/* 3D double-precision, complex evaulation functions        */
/************************************************************/
void
eval_multi_UBspline_3d_z (multi_UBspline_3d_z *spline,
                          double x, double y, double z,
                          complex_double* restrict vals)
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
  complex_double* restrict coefs = spline->coefs;
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
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      vector double abc0 = vec_splats(a[i]*b[j]*c[0]);
      vector double abc1 = vec_splats(a[i]*b[j]*c[1]);
      vector double abc2 = vec_splats(a[i]*b[j]*c[2]);
      vector double abc3 = vec_splats(a[i]*b[j]*c[3]);
      // double abc0 = a[i]*b[j]*c[0];
      // double abc1 = a[i]*b[j]*c[1];
      // double abc2 = a[i]*b[j]*c[2];
      // double abc3 = a[i]*b[j]*c[3];
//       complex_double* restrict coefs0 = spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+0)*zs);
//       complex_double* restrict coefs1 = spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+1)*zs);
//       complex_double* restrict coefs2 = spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+2)*zs);
//       complex_double* restrict coefs3 = spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+3)*zs);
      vector double* restrict coefs0 = (vector double*)(spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+0)*zs));
      vector double* restrict coefs1 = (vector double*)(spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+1)*zs));
      vector double* restrict coefs2 = (vector double*)(spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+2)*zs));
      vector double* restrict coefs3 = (vector double*)(spline->coefs + ((ix+i)*xs +(iy+j)*ys +(iz+3)*zs));
      for (int n=0; n<(spline->num_splines-1); n+=2)
      {
        vector double v0 = *((vector double*)&(vals[n]));
        vector double v1 = *((vector double*)&(vals[n+1]));
        v0 = vec_madd (coefs0[n  ], abc0, v0);
        v0 = vec_madd (coefs1[n  ], abc1, v0);
        v0 = vec_madd (coefs2[n  ], abc2, v0);
        v0 = vec_madd (coefs3[n  ], abc3, v0);
        v1 = vec_madd (coefs0[n+1], abc0, v1);
        v1 = vec_madd (coefs1[n+1], abc1, v1);
        v1 = vec_madd (coefs2[n+1], abc2, v1);
        v1 = vec_madd (coefs3[n+1], abc3, v1);
        *((vector double*)&(vals[n]))   = v0;
        *((vector double*)&(vals[n+1])) = v1;
      }
      if (spline->num_splines & 1)
      {
        int n = spline->num_splines-1;
        vector double v0 = *((vector double*)&(vals[n]));
        v0 = vec_madd (coefs0[n  ], abc0, v0);
        v0 = vec_madd (coefs1[n  ], abc1, v0);
        v0 = vec_madd (coefs2[n  ], abc2, v0);
        v0 = vec_madd (coefs3[n  ], abc3, v0);
        *((vector double*)&(vals[n]))   = v0;
      }
    }
}


void
eval_multi_UBspline_3d_z_vg (multi_UBspline_3d_z *spline,
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
        complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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
eval_multi_UBspline_3d_z_vgl (multi_UBspline_3d_z *spline,
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
        complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
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
eval_multi_UBspline_3d_z_vgh (multi_UBspline_3d_z *spline,
                              double x, double y, double z,
                              complex_double* /*restrict*/ vals,
                              complex_double* /*restrict*/ grads,
                              complex_double* /*restrict*/ hess)
{
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }
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
  complex_double* restrict coefs = spline->coefs;
  double dxInv = spline->x_grid.delta_inv;
  double dyInv = spline->y_grid.delta_inv;
  double dzInv = spline->z_grid.delta_inv;
  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = dxInv*(dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = dxInv*(dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = dxInv*(dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = dxInv*(dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = dxInv*dxInv*(d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = dxInv*dxInv*(d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = dxInv*dxInv*(d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = dxInv*dxInv*(d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);
  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = dyInv*(dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = dyInv*(dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = dyInv*(dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = dyInv*(dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = dyInv*dyInv*(d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = dyInv*dyInv*(d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = dyInv*dyInv*(d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = dyInv*dyInv*(d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);
  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = dzInv*(dAd[ 0]*tpz[0] + dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = dzInv*(dAd[ 4]*tpz[0] + dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = dzInv*(dAd[ 8]*tpz[0] + dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = dzInv*(dAd[12]*tpz[0] + dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = dzInv*dzInv*(d2Ad[ 0]*tpz[0] + d2Ad[ 1]*tpz[1] + d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = dzInv*dzInv*(d2Ad[ 4]*tpz[0] + d2Ad[ 5]*tpz[1] + d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = dzInv*dzInv*(d2Ad[ 8]*tpz[0] + d2Ad[ 9]*tpz[1] + d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = dzInv*dzInv*(d2Ad[12]*tpz[0] + d2Ad[13]*tpz[1] + d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);
  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      //      	vector double abc_v[40];
      vector double abc_0 = vec_splats(  a[i] *  b[j] *  c[0]);
      vector double abc_1 = vec_splats(  a[i] *  b[j] *  c[1]);
      vector double abc_2 = vec_splats(  a[i] *  b[j] *  c[2]);
      vector double abc_3 = vec_splats(  a[i] *  b[j] *  c[3]);
      vector double abc_4 = vec_splats( da[i] *  b[j] *  c[0]);
      vector double abc_5 = vec_splats( da[i] *  b[j] *  c[1]);
      vector double abc_6 = vec_splats( da[i] *  b[j] *  c[2]);
      vector double abc_7 = vec_splats( da[i] *  b[j] *  c[3]);
      vector double abc_8 = vec_splats(  a[i] * db[j]  *  c[0]);
      vector double abc_9 = vec_splats(  a[i] * db[j]  *  c[1]);
      vector double abc_10 = vec_splats(  a[i] * db[j] *  c[2]);
      vector double abc_11 = vec_splats(  a[i] * db[j] *  c[3]);
      vector double abc_12 = vec_splats(  a[i] *  b[j] * dc[0]);
      vector double abc_13 = vec_splats(  a[i] *  b[j] * dc[1]);
      vector double abc_14 = vec_splats(  a[i] *  b[j] * dc[2]);
      vector double abc_15 = vec_splats(  a[i] *  b[j] * dc[3]);
      vector double abc_16 = vec_splats(d2a[i] *  b[j] *  c[0]);
      vector double abc_17 = vec_splats(d2a[i] *  b[j] *  c[1]);
      vector double abc_18 = vec_splats(d2a[i] *  b[j] *  c[2]);
      vector double abc_19 = vec_splats(d2a[i] *  b[j] *  c[3]);
      vector double abc_20 = vec_splats( da[i] * db[j] *  c[0]);
      vector double abc_21 = vec_splats( da[i] * db[j] *  c[1]);
      vector double abc_22 = vec_splats( da[i] * db[j] *  c[2]);
      vector double abc_23 = vec_splats( da[i] * db[j] *  c[3]);
      vector double abc_24 = vec_splats( da[i] *  b[j] * dc[0]);
      vector double abc_25 = vec_splats( da[i] *  b[j] * dc[1]);
      vector double abc_26 = vec_splats( da[i] *  b[j] * dc[2]);
      vector double abc_27 = vec_splats( da[i] *  b[j] * dc[3]);
      vector double abc_28 = vec_splats(  a[i] *d2b[j] *  c[0]);
      vector double abc_29 = vec_splats(  a[i] *d2b[j] *  c[1]);
      vector double abc_30 = vec_splats(  a[i] *d2b[j] *  c[2]);
      vector double abc_31 = vec_splats(  a[i] *d2b[j] *  c[3]);
      vector double abc_32 = vec_splats(  a[i] * db[j] * dc[0]);
      vector double abc_33 = vec_splats(  a[i] * db[j] * dc[1]);
      vector double abc_34 = vec_splats(  a[i] * db[j] * dc[2]);
      vector double abc_35 = vec_splats(  a[i] * db[j] * dc[3]);
      vector double abc_36 = vec_splats(  a[i] *  b[j] *d2c[0]);
      vector double abc_37 = vec_splats(  a[i] *  b[j] *d2c[1]);
      vector double abc_38 = vec_splats(  a[i] *  b[j] *d2c[2]);
      vector double abc_39 = vec_splats(  a[i] *  b[j] *d2c[3]);
      complex_double* restrict coefs =
        (complex_double* restrict)spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz)*zs);
      for (int n=0; n<spline->num_splines; n++)
      {
        vector double c0 = *(vector double*)&(coefs[n     ]);//vec_xld2 (n     ,coefs);
        vector double c1 = *(vector double*)&(coefs[n+  zs]);//vec_xld2 (n+  zs,coefs);
        vector double c2 = *(vector double*)&(coefs[n+2*zs]);//vec_xld2 (n+2*zs,coefs);
        vector double c3 = *(vector double*)&(coefs[n+3*zs]);//vec_xld2 (n+3*zs,coefs);
        vector double v = *(vector double*)&(vals[n]);
        v           = vec_madd (c0, abc_0, v);
        v           = vec_madd (c1, abc_1, v);
        v           = vec_madd (c2, abc_2, v);
        v           = vec_madd (c3, abc_3, v);
        *(vector double*)&(vals[n]) = v;
        vector double g0 = *(vector double*)&(grads[3*n+0]);
        vector double g1 = *(vector double*)&(grads[3*n+1]);
        vector double g2 = *(vector double*)&(grads[3*n+2]);
        g0 = vec_madd (c0, abc_4 , g0);
        g0 = vec_madd (c1, abc_5 , g0);
        g0 = vec_madd (c2, abc_6 , g0);
        g0 = vec_madd (c3, abc_7 , g0);
        g1 = vec_madd (c0, abc_8 , g1);
        g1 = vec_madd (c1, abc_9 , g1);
        g1 = vec_madd (c2, abc_10, g1);
        g1 = vec_madd (c3, abc_11, g1);
        g2 = vec_madd (c0, abc_12, g2);
        g2 = vec_madd (c1, abc_13, g2);
        g2 = vec_madd (c2, abc_14, g2);
        g2 = vec_madd (c3, abc_15, g2);
        *(vector double*)&(grads[3*n+0]) = g0;
        *(vector double*)&(grads[3*n+1]) = g1;
        *(vector double*)&(grads[3*n+2]) = g2;
        vector double h0 = *(vector double*)&(hess[9*n+0]);
        vector double h1 = *(vector double*)&(hess[9*n+1]);
        vector double h2 = *(vector double*)&(hess[9*n+2]);
        vector double h4 = *(vector double*)&(hess[9*n+4]);
        vector double h5 = *(vector double*)&(hess[9*n+5]);
        vector double h8 = *(vector double*)&(hess[9*n+8]);
        h0 = vec_madd (c0, abc_16, h0);
        h0 = vec_madd (c1, abc_17, h0);
        h0 = vec_madd (c2, abc_18, h0);
        h0 = vec_madd (c3, abc_19, h0);
        h1 = vec_madd (c0, abc_20, h1);
        h1 = vec_madd (c1, abc_21, h1);
        h1 = vec_madd (c2, abc_22, h1);
        h1 = vec_madd (c3, abc_23, h1);
        h2 = vec_madd (c0, abc_24, h2);
        h2 = vec_madd (c1, abc_25, h2);
        h2 = vec_madd (c2, abc_26, h2);
        h2 = vec_madd (c3, abc_27, h2);
        h4 = vec_madd (c0, abc_28, h4);
        h4 = vec_madd (c1, abc_29, h4);
        h4 = vec_madd (c2, abc_30, h4);
        h4 = vec_madd (c3, abc_31, h4);
        h5 = vec_madd (c0, abc_32, h5);
        h5 = vec_madd (c1, abc_33, h5);
        h5 = vec_madd (c2, abc_34, h5);
        h5 = vec_madd (c3, abc_35, h5);
        h8 = vec_madd (c0, abc_36, h8);
        h8 = vec_madd (c1, abc_37, h8);
        h8 = vec_madd (c2, abc_38, h8);
        h8 = vec_madd (c3, abc_39, h8);
        *(vector double*)&(hess[9*n+0]) = h0;
        *(vector double*)&(hess[9*n+1]) = h1;
        *(vector double*)&(hess[9*n+2]) = h2;
        *(vector double*)&(hess[9*n+4]) = h4;
        *(vector double*)&(hess[9*n+5]) = h5;
        *(vector double*)&(hess[9*n+8]) = h8;
      }
    }
  for (int n=0; n<spline->num_splines; n++)
  {
    // Copy hessian elements into lower half of 3x3 matrix
    hess[9*n+3] = hess[9*n+1];
    hess[9*n+6] = hess[9*n+2];
    hess[9*n+7] = hess[9*n+5];
  }
}

// void
// eval_multi_UBspline_3d_z_vgh (multi_UBspline_3d_z *spline,
// 			      double x, double y, double z,
// 			      complex_double* restrict vals,
// 			      complex_double* restrict grads,
// 			      complex_double* restrict hess)
// {
//   x -= spline->x_grid.start;
//   y -= spline->y_grid.start;
//   z -= spline->z_grid.start;
//   double ux = x*spline->x_grid.delta_inv;
//   double uy = y*spline->y_grid.delta_inv;
//   double uz = z*spline->z_grid.delta_inv;
//   double ipartx, iparty, ipartz, tx, ty, tz;
//   tx = modf (ux, &ipartx);  int ix = (int) ipartx;
//   ty = modf (uy, &iparty);  int iy = (int) iparty;
//   tz = modf (uz, &ipartz);  int iz = (int) ipartz;

//   double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4],
//     da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
//   tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
//   tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
//   tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
//   complex_double* restrict coefs = spline->coefs;

//   a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
//   a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
//   a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
//   a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
//   da[0] = (dAd[ 0]*tpx[0] + dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
//   da[1] = (dAd[ 4]*tpx[0] + dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
//   da[2] = (dAd[ 8]*tpx[0] + dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
//   da[3] = (dAd[12]*tpx[0] + dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);
//   d2a[0] = (d2Ad[ 0]*tpx[0] + d2Ad[ 1]*tpx[1] + d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
//   d2a[1] = (d2Ad[ 4]*tpx[0] + d2Ad[ 5]*tpx[1] + d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
//   d2a[2] = (d2Ad[ 8]*tpx[0] + d2Ad[ 9]*tpx[1] + d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
//   d2a[3] = (d2Ad[12]*tpx[0] + d2Ad[13]*tpx[1] + d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

//   b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
//   b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
//   b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
//   b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
//   db[0] = (dAd[ 0]*tpy[0] + dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
//   db[1] = (dAd[ 4]*tpy[0] + dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
//   db[2] = (dAd[ 8]*tpy[0] + dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
//   db[3] = (dAd[12]*tpy[0] + dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
//   d2b[0] = (d2Ad[ 0]*tpy[0] + d2Ad[ 1]*tpy[1] + d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
//   d2b[1] = (d2Ad[ 4]*tpy[0] + d2Ad[ 5]*tpy[1] + d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
//   d2b[2] = (d2Ad[ 8]*tpy[0] + d2Ad[ 9]*tpy[1] + d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
//   d2b[3] = (d2Ad[12]*tpy[0] + d2Ad[13]*tpy[1] + d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

//   c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
//   c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
//   c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
//   c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
//   dc[0] = (dAd[ 0]*tpz[0] + dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
//   dc[1] = (dAd[ 4]*tpz[0] + dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
//   dc[2] = (dAd[ 8]*tpz[0] + dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
//   dc[3] = (dAd[12]*tpz[0] + dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
//   d2c[0] = (d2Ad[ 0]*tpz[0] + d2Ad[ 1]*tpz[1] + d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
//   d2c[1] = (d2Ad[ 4]*tpz[0] + d2Ad[ 5]*tpz[1] + d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
//   d2c[2] = (d2Ad[ 8]*tpz[0] + d2Ad[ 9]*tpz[1] + d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
//   d2c[3] = (d2Ad[12]*tpz[0] + d2Ad[13]*tpz[1] + d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

//   intptr_t xs = spline->x_stride;
//   intptr_t ys = spline->y_stride;
//   intptr_t zs = spline->z_stride;

//   for (int n=0; n<spline->num_splines; n++) {
//     vals[n] = 0.0;
//     grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
//     for (int i=0; i<9; i++)
//       hess[9*n+i] = 0.0;
//   }

//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       //      	vector double abc_v[40];
// 	vector double abc_0 = vec_splats(  a[i] *  b[j] *  c[0]);
// 	vector double abc_1 = vec_splats(  a[i] *  b[j] *  c[1]);
// 	vector double abc_2 = vec_splats(  a[i] *  b[j] *  c[2]);
// 	vector double abc_3 = vec_splats(  a[i] *  b[j] *  c[3]);

// 	vector double abc_4 = vec_splats( da[i] *  b[j] *  c[0]);
// 	vector double abc_5 = vec_splats( da[i] *  b[j] *  c[1]);
// 	vector double abc_6 = vec_splats( da[i] *  b[j] *  c[2]);
// 	vector double abc_7 = vec_splats( da[i] *  b[j] *  c[3]);

// 	vector double abc_8 = vec_splats(  a[i] * db[j]  *  c[0]);
// 	vector double abc_9 = vec_splats(  a[i] * db[j]  *  c[1]);
// 	vector double abc_10 = vec_splats(  a[i] * db[j] *  c[2]);
// 	vector double abc_11 = vec_splats(  a[i] * db[j] *  c[3]);

// 	vector double abc_12 = vec_splats(  a[i] *  b[j] * dc[0]);
// 	vector double abc_13 = vec_splats(  a[i] *  b[j] * dc[1]);
// 	vector double abc_14 = vec_splats(  a[i] *  b[j] * dc[2]);
// 	vector double abc_15 = vec_splats(  a[i] *  b[j] * dc[3]);

// 	vector double abc_16 = vec_splats(d2a[i] *  b[j] *  c[0]);
// 	vector double abc_17 = vec_splats(d2a[i] *  b[j] *  c[1]);
// 	vector double abc_18 = vec_splats(d2a[i] *  b[j] *  c[2]);
// 	vector double abc_19 = vec_splats(d2a[i] *  b[j] *  c[3]);


// 	double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz)*zs);
// 	for (int n=0; n<spline->num_splines; n++) {
// 	  vector double c0 = *(vector double*)&(coefs[n     ]);//vec_xld2 (n     ,coefs);
// 	  vector double c1 = *(vector double*)&(coefs[n+  zs]);//vec_xld2 (n+  zs,coefs);
// 	  vector double c2 = *(vector double*)&(coefs[n+2*zs]);//vec_xld2 (n+2*zs,coefs);
// 	  vector double c3 = *(vector double*)&(coefs[n+3*zs]);//vec_xld2 (n+3*zs,coefs);

// 	  vector double v = *(vector double*)&(vals[n]);
// 	  v           = vec_madd (c0, abc_0, v);
// 	  v           = vec_madd (c1, abc_1, v);
// 	  v           = vec_madd (c2, abc_2, v);
// 	  v           = vec_madd (c3, abc_3, v);
// 	  *(vector double*)&(vals[n]) = v;
// 	  vector double g0 = *(vector double*)&(grads[3*n+0]);
// 	  vector double g1 = *(vector double*)&(grads[3*n+1]);
// 	  vector double g2 = *(vector double*)&(grads[3*n+2]);
// 	  g0 = vec_madd (c0, abc_4, g0);
// 	  g0 = vec_madd (c1, abc_5, g0);
// 	  g0 = vec_madd (c2, abc_6, g0);
// 	  g0 = vec_madd (c3, abc_7, g0);
// 	  g1 = vec_madd (c0, abc_8, g1);
// 	  g1 = vec_madd (c1, abc_9, g1);
// 	  g1 = vec_madd (c2, abc_10, g1);
// 	  g1 = vec_madd (c3, abc_11, g1);
// 	  g2 = vec_madd (c0, abc_12, g2);
// 	  g2 = vec_madd (c1, abc_13, g2);
// 	  g2 = vec_madd (c2, abc_14, g2);
// 	  g2 = vec_madd (c3, abc_15, g2);
// 	  *(vector double*)&(grads[3*n+0]) = g0;
// 	  *(vector double*)&(grads[3*n+2]) = g1;
// 	  *(vector double*)&(grads[3*n+4]) = g2;
// 	  vector double h0 = *(vector double*)&(hess[9*n+0]);
// 	  h0 = vec_madd (c0, abc_16, h0);
// 	  h0 = vec_madd (c1, abc_17, h0);
// 	  h0 = vec_madd (c2, abc_18, h0);
// 	  h0 = vec_madd (c3, abc_19, h0);
// 	  *(vector double*)&(hess[9*n+0]) = h0;
// 	}
//     }

//   for (int i=0; i<4; i++)
//     for (int j=0; j<4; j++) {
//       vector double abc_20 = vec_splats( da[i] * db[j] *  c[0]);
//       vector double abc_21 = vec_splats( da[i] * db[j] *  c[1]);
//       vector double abc_22 = vec_splats( da[i] * db[j] *  c[2]);
//       vector double abc_23 = vec_splats( da[i] * db[j] *  c[3]);

//       vector double abc_24 = vec_splats( da[i] *  b[j] * dc[0]);
//       vector double abc_25 = vec_splats( da[i] *  b[j] * dc[1]);
//       vector double abc_26 = vec_splats( da[i] *  b[j] * dc[2]);
//       vector double abc_27 = vec_splats( da[i] *  b[j] * dc[3]);

//       vector double abc_28 = vec_splats(  a[i] *d2b[j] *  c[0]);
//       vector double abc_29 = vec_splats(  a[i] *d2b[j] *  c[1]);
//       vector double abc_30 = vec_splats(  a[i] *d2b[j] *  c[2]);
//       vector double abc_31 = vec_splats(  a[i] *d2b[j] *  c[3]);

//       vector double abc_32 = vec_splats(  a[i] * db[j] * dc[0]);
//       vector double abc_33 = vec_splats(  a[i] * db[j] * dc[1]);
//       vector double abc_34 = vec_splats(  a[i] * db[j] * dc[2]);
//       vector double abc_35 = vec_splats(  a[i] * db[j] * dc[3]);

//       vector double abc_36 = vec_splats(  a[i] *  b[j] *d2c[0]);
//       vector double abc_37 = vec_splats(  a[i] *  b[j] *d2c[1]);
//       vector double abc_38 = vec_splats(  a[i] *  b[j] *d2c[2]);
//       vector double abc_39 = vec_splats(  a[i] *  b[j] *d2c[3]);

//       double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz)*zs);
//       for (int n=0; n<spline->num_splines; n++) {
// 	vector double c0 = *(vector double*)&(coefs[n     ]);//vec_xld2 (n     ,coefs);
// 	vector double c1 = *(vector double*)&(coefs[n+  zs]);//vec_xld2 (n+  zs,coefs);
// 	vector double c2 = *(vector double*)&(coefs[n+2*zs]);//vec_xld2 (n+2*zs,coefs);
// 	vector double c3 = *(vector double*)&(coefs[n+3*zs]);//vec_xld2 (n+3*zs,coefs);

// 	vector double h1 = *(vector double*)&(hess[9*n+1]);
// 	vector double h2 = *(vector double*)&(hess[9*n+2]);
// 	vector double h4 = *(vector double*)&(hess[9*n+4]);
// 	vector double h5 = *(vector double*)&(hess[9*n+5]);
// 	vector double h8 = *(vector double*)&(hess[9*n+8]);

// 	h1 = vec_madd (c0, abc_20, h1);
// 	h1 = vec_madd (c1, abc_21, h1);
// 	h1 = vec_madd (c2, abc_22, h1);
// 	h1 = vec_madd (c3, abc_23, h1);
// 	h2 = vec_madd (c0, abc_24, h2);
// 	h2 = vec_madd (c1, abc_25, h2);
// 	h2 = vec_madd (c2, abc_26, h2);
// 	h2 = vec_madd (c3, abc_27, h2);
// 	h4 = vec_madd (c0, abc_28, h4);
// 	h4 = vec_madd (c1, abc_29, h4);
// 	h4 = vec_madd (c2, abc_30, h4);
// 	h4 = vec_madd (c3, abc_31, h4);
// 	h5 = vec_madd (c0, abc_32, h5);
// 	h5 = vec_madd (c1, abc_33, h5);
// 	h5 = vec_madd (c2, abc_34, h5);
// 	h5 = vec_madd (c3, abc_35, h5);
// 	h8 = vec_madd (c0, abc_36, h8);
// 	h8 = vec_madd (c1, abc_37, h8);
// 	h8 = vec_madd (c2, abc_38, h8);
// 	h8 = vec_madd (c3, abc_39, h8);

// 	*(vector double*)&(hess[9*n+1]) = h1;
// 	*(vector double*)&(hess[9*n+2]) = h2;
// 	*(vector double*)&(hess[9*n+4]) = h4;
// 	*(vector double*)&(hess[9*n+5]) = h5;
// 	*(vector double*)&(hess[9*n+8]) = h8;

//       }
//     }

//   double dxInv = spline->x_grid.delta_inv;
//   double dyInv = spline->y_grid.delta_inv;
//   double dzInv = spline->z_grid.delta_inv;
//   for (int n=0; n<spline->num_splines; n++) {
//     grads[3*n+0] *= dxInv;
//     grads[3*n+1] *= dyInv;
//     grads[3*n+2] *= dzInv;
//     hess[9*n+0] *= dxInv*dxInv;
//     hess[9*n+4] *= dyInv*dyInv;
//     hess[9*n+8] *= dzInv*dzInv;
//     hess[9*n+1] *= dxInv*dyInv;
//     hess[9*n+2] *= dxInv*dzInv;
//     hess[9*n+5] *= dyInv*dzInv;
//     // Copy hessian elements into lower half of 3x3 matrix
//     hess[9*n+3] = hess[9*n+1];
//     hess[9*n+6] = hess[9*n+2];
//     hess[9*n+7] = hess[9*n+5];
//   }
// }




#endif
