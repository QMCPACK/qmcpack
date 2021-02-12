/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "bspline_base.h"
#include "multi_nubspline_structs.h"
#include "multi_nubspline_eval_z.h"

/************************************************************/
/* 1D double-precision, complex evaluation functions        */
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
