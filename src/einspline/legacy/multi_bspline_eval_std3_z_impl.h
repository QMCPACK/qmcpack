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

#ifndef MULTI_BSPLINE_EVAL_STD3_Z_IMPL_H
#define MULTI_BSPLINE_EVAL_STD3_Z_IMPL_H

#include <math.h>
#include <stdio.h>
#include <einspline/bspline_base.h>
#include <einspline/multi_bspline_structs.h>

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
eval_multi_UBspline_1d_z_vg (const multi_UBspline_1d_z *spline,
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
eval_multi_UBspline_1d_z_vgl (const multi_UBspline_1d_z *spline,
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
eval_multi_UBspline_3d_z (const multi_UBspline_3d_z *spline,
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
      for (int k=0; k<4; k++)
      {
        double prefactor = a[i]*b[j]*c[k];
        complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<spline->num_splines; n++)
          vals[n] += prefactor*coefs[n];
      }
}



#if 1
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
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
  /*
  	      double pre0 =   a[i] *   b[j];
  	      double pre1 =  da[i] *   b[j];
  	      double pre2 = d2a[i] *   b[j];
                double pre3 =   a[i] *  db[j];
                double pre4 =  da[i] *  db[j];
  	      double pre5 =   a[i] * d2b[j];

  	      for (int k=0; k<4; k++) {
  		  complex_double coef = coefs[i*xs + j*ys + k*zs];
                    complex_double   ck = c[k] * coef;
                    complex_double  dck = dc[k] * coef;
  		  val   +=   pre0 *  ck;
  		  grad0 +=   pre1 *  ck;
  		  grad1 +=   pre3 *  ck;
  		  grad2 +=   pre0 * dck;
  		  hess0 +=   pre2 *  ck;
  		  hess1 +=   pre4 *  ck;
  		  hess2 +=   pre1 * dck;
  		  hess4 +=   pre5 *  ck;
  		  hess5 +=   pre3 * dck;
  		  hess8 +=   pre0 *d2c[k] * coef;
  	      }
  */
  complex_double* restrict coefs __attribute__((__aligned__(16))) = spline->coefs + ix*xs + iy*ys + iz*zs;
  int num_splines = spline->num_splines;
  for (int n=0; n<num_splines; n++, coefs++)
  {
    complex_double val = 0, grad0 = 0, grad1 = 0, grad2 = 0;
    complex_double hess0 = 0, hess1 = 0, hess2 = 0, hess4 = 0, hess5 = 0, hess8 = 0;
    for (int i=0; i<4; i++)
    {
#if 0
#pragma unroll(4)
      for (int j=0; j<4; j++)
      {
        complex_double* restrict coefp __attribute__((__aligned__(16))) = coefs + i*xs + j*ys;
        double pre0 =   a[i] *   b[j];
        double pre1 =  da[i] *   b[j];
        double pre2 = d2a[i] *   b[j];
        double pre3 =   a[i] *  db[j];
        complex_double coef0 = coefp[0];
        complex_double coef1 = coefp[zs];
        complex_double coef2 = coefp[2*zs];
        complex_double coef3 = coefp[3*zs];
        complex_double sum0 = c[0] * coef0 + c[1] * coef1 + c[2] * coef2 + c[3] * coef3;
        complex_double sum1 = dc[0] * coef0 + dc[1] * coef1 + dc[2] * coef2 + dc[3] * coef3;
        hess8 +=   pre0 * (d2c[0] * coef0 +d2c[1] * coef1 + d2c[2] * coef2  + d2c[3] * coef3);
        hess1 +=   (da[i] *  db[j]) * sum0;
        hess4 +=   (a[i] * d2b[j]) * sum0;
        val   +=   pre0 * sum0;
        grad0 +=   pre1 * sum0;
        grad1 +=   pre3 * sum0;
        grad2 +=   pre0 * sum1;
        hess0 +=   pre2 * sum0;
        hess2 +=   pre1 * sum1;
        hess5 +=   pre3 * sum1;
      }
#else
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
      hess0 +=   pre21 * sum01 + pre2 * sum0 + pre22 * sum02 + pre23 * sum03;
#endif
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
  for (int n=0; n<num_splines; n++)
  {
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
#else
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
        double abc = a[i]*b[j]*c[k];
        double dabc[3], d2abc[6];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];
        d2abc[0] = d2a[i]*  b[j]*  c[k];
        d2abc[1] =  da[i]* db[j]*  c[k];
        d2abc[2] =  da[i]*  b[j]* dc[k];
        d2abc[3] =   a[i]*d2b[j]*  c[k];
        d2abc[4] =   a[i]* db[j]* dc[k];
        d2abc[5] =   a[i]*  b[j]*d2c[k];
        complex_double* restrict coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (int n=0; n<spline->num_splines; n++)
        {
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
  tx = modf (ux, &ipartx);
  int ix = (int) ipartx;
  ty = modf (uy, &iparty);
  int iy = (int) iparty;
  tz = modf (uz, &ipartz);
  int iz = (int) ipartz;
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4],
         da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4],
         d3a[4], d3b[4], d3c[4];
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
  for (int n=0; n<spline->num_splines; n++)
  {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (int i=0; i<9; i++)
      hess[9*n+i] = 0.0;
    for (int i=0; i<27; i++)
      gradhess[27*n+i] = 0.0;
  }
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++)
      {
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
        for (int n=0; n<spline->num_splines; n++)
        {
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
