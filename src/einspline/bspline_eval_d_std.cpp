/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include "bspline_base.h"
#include "bspline_structs.h"
#include "bspline_eval_d.h"

extern const double* restrict Ad;
extern const double* restrict dAd;
extern const double* restrict d2Ad;

/** get the i offset and the t remainder
 *  for particular grid and boundary condition
 *  if debug check bounds
 */
inline void coeff_offset_UBspline_d(const Ugrid* const restrict grid,
                                    const BCtype_d* const restrict BC,
                                    double x,
                                    int& i,
                                    double& t)
{
  x -= grid->start;
  double u = x * grid->delta_inv;
  double ipart;
  t = modf(u, &ipart);
  i = ipart;
#ifndef NDEBUG
  assert(ipart >= 0);
  if (BC->lCode == NATURAL)
  {
    assert(ipart <= grid->num - 2);
  }
  else
  {
    assert(ipart <= grid->num - 1);
  }
#endif
}


/************************************************************/
/* 1D double-precision, real evaluation functions           */
/************************************************************/


/** Value only */
void eval_UBspline_1d_d(UBspline_1d_d* restrict spline, double x, double* restrict val)
{
  int i;
  double t;
  coeff_offset_UBspline_d(&(spline->x_grid), &(spline->xBC), x, i, t);
  double tp[4];
  tp[0]                  = t * t * t;
  tp[1]                  = t * t;
  tp[2]                  = t;
  tp[3]                  = 1.0;
  double* restrict coefs = spline->coefs;
  *val                   = coefs[i + 0] * (Ad[0] * tp[0] + Ad[1] * tp[1] + Ad[2] * tp[2] + Ad[3] * tp[3]);
  *val += coefs[i + 1] * (Ad[4] * tp[0] + Ad[5] * tp[1] + Ad[6] * tp[2] + Ad[7] * tp[3]);
  *val += coefs[i + 2] * (Ad[8] * tp[0] + Ad[9] * tp[1] + Ad[10] * tp[2] + Ad[11] * tp[3]);
  *val += coefs[i + 3] * (Ad[12] * tp[0] + Ad[13] * tp[1] + Ad[14] * tp[2] + Ad[15] * tp[3]);
}

/************************************************************/
/* 3D double-precision, real evaluation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_d(UBspline_3d_d* restrict spline, double x, double y, double z, double* restrict val)
{
  int ix;
  int iy;
  int iz;
  double tx;
  double ty;
  double tz;
  coeff_offset_UBspline_d(&(spline->x_grid), &(spline->xBC), x, ix, tx);
  coeff_offset_UBspline_d(&(spline->y_grid), &(spline->yBC), y, iy, ty);
  coeff_offset_UBspline_d(&(spline->z_grid), &(spline->zBC), z, iz, tz);
  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0]                 = tx * tx * tx;
  tpx[1]                 = tx * tx;
  tpx[2]                 = tx;
  tpx[3]                 = 1.0;
  tpy[0]                 = ty * ty * ty;
  tpy[1]                 = ty * ty;
  tpy[2]                 = ty;
  tpy[3]                 = 1.0;
  tpz[0]                 = tz * tz * tz;
  tpz[1]                 = tz * tz;
  tpz[2]                 = tz;
  tpz[3]                 = 1.0;
  double* restrict coefs = spline->coefs;
  a[0]                   = (Ad[0] * tpx[0] + Ad[1] * tpx[1] + Ad[2] * tpx[2] + Ad[3] * tpx[3]);
  a[1]                   = (Ad[4] * tpx[0] + Ad[5] * tpx[1] + Ad[6] * tpx[2] + Ad[7] * tpx[3]);
  a[2]                   = (Ad[8] * tpx[0] + Ad[9] * tpx[1] + Ad[10] * tpx[2] + Ad[11] * tpx[3]);
  a[3]                   = (Ad[12] * tpx[0] + Ad[13] * tpx[1] + Ad[14] * tpx[2] + Ad[15] * tpx[3]);
  b[0]                   = (Ad[0] * tpy[0] + Ad[1] * tpy[1] + Ad[2] * tpy[2] + Ad[3] * tpy[3]);
  b[1]                   = (Ad[4] * tpy[0] + Ad[5] * tpy[1] + Ad[6] * tpy[2] + Ad[7] * tpy[3]);
  b[2]                   = (Ad[8] * tpy[0] + Ad[9] * tpy[1] + Ad[10] * tpy[2] + Ad[11] * tpy[3]);
  b[3]                   = (Ad[12] * tpy[0] + Ad[13] * tpy[1] + Ad[14] * tpy[2] + Ad[15] * tpy[3]);
  c[0]                   = (Ad[0] * tpz[0] + Ad[1] * tpz[1] + Ad[2] * tpz[2] + Ad[3] * tpz[3]);
  c[1]                   = (Ad[4] * tpz[0] + Ad[5] * tpz[1] + Ad[6] * tpz[2] + Ad[7] * tpz[3]);
  c[2]                   = (Ad[8] * tpz[0] + Ad[9] * tpz[1] + Ad[10] * tpz[2] + Ad[11] * tpz[3]);
  c[3]                   = (Ad[12] * tpz[0] + Ad[13] * tpz[1] + Ad[14] * tpz[2] + Ad[15] * tpz[3]);
  int xs                 = spline->x_stride;
  int ys                 = spline->y_stride;
#define P(i, j, k) coefs[(ix + (i)) * xs + (iy + (j)) * ys + (iz + (k))]
  *val = (a[0] *
              (b[0] * (P(0, 0, 0) * c[0] + P(0, 0, 1) * c[1] + P(0, 0, 2) * c[2] + P(0, 0, 3) * c[3]) +
               b[1] * (P(0, 1, 0) * c[0] + P(0, 1, 1) * c[1] + P(0, 1, 2) * c[2] + P(0, 1, 3) * c[3]) +
               b[2] * (P(0, 2, 0) * c[0] + P(0, 2, 1) * c[1] + P(0, 2, 2) * c[2] + P(0, 2, 3) * c[3]) +
               b[3] * (P(0, 3, 0) * c[0] + P(0, 3, 1) * c[1] + P(0, 3, 2) * c[2] + P(0, 3, 3) * c[3])) +
          a[1] *
              (b[0] * (P(1, 0, 0) * c[0] + P(1, 0, 1) * c[1] + P(1, 0, 2) * c[2] + P(1, 0, 3) * c[3]) +
               b[1] * (P(1, 1, 0) * c[0] + P(1, 1, 1) * c[1] + P(1, 1, 2) * c[2] + P(1, 1, 3) * c[3]) +
               b[2] * (P(1, 2, 0) * c[0] + P(1, 2, 1) * c[1] + P(1, 2, 2) * c[2] + P(1, 2, 3) * c[3]) +
               b[3] * (P(1, 3, 0) * c[0] + P(1, 3, 1) * c[1] + P(1, 3, 2) * c[2] + P(1, 3, 3) * c[3])) +
          a[2] *
              (b[0] * (P(2, 0, 0) * c[0] + P(2, 0, 1) * c[1] + P(2, 0, 2) * c[2] + P(2, 0, 3) * c[3]) +
               b[1] * (P(2, 1, 0) * c[0] + P(2, 1, 1) * c[1] + P(2, 1, 2) * c[2] + P(2, 1, 3) * c[3]) +
               b[2] * (P(2, 2, 0) * c[0] + P(2, 2, 1) * c[1] + P(2, 2, 2) * c[2] + P(2, 2, 3) * c[3]) +
               b[3] * (P(2, 3, 0) * c[0] + P(2, 3, 1) * c[1] + P(2, 3, 2) * c[2] + P(2, 3, 3) * c[3])) +
          a[3] *
              (b[0] * (P(3, 0, 0) * c[0] + P(3, 0, 1) * c[1] + P(3, 0, 2) * c[2] + P(3, 0, 3) * c[3]) +
               b[1] * (P(3, 1, 0) * c[0] + P(3, 1, 1) * c[1] + P(3, 1, 2) * c[2] + P(3, 1, 3) * c[3]) +
               b[2] * (P(3, 2, 0) * c[0] + P(3, 2, 1) * c[1] + P(3, 2, 2) * c[2] + P(3, 2, 3) * c[3]) +
               b[3] * (P(3, 3, 0) * c[0] + P(3, 3, 1) * c[1] + P(3, 3, 2) * c[2] + P(3, 3, 3) * c[3])));
#undef P
}
