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

#ifndef BSPLINE_EVAL_SSE_S_H
#define BSPLINE_EVAL_SSE_S_H

#include <stdio.h>
#include <math.h>
#include <ppc_intrinsics.h>

extern vector float  A0,   A1,   A2,   A3;
extern vector float  dA0,  dA1,  dA2,  dA3;
extern vector float d2A0, d2A1, d2A2, d2A3;

extern const float* restrict   Af;
extern const float* restrict  dAf;
extern const float* restrict d2Af;

inline vector float
MakeVec (double a, double b, double c, double d)
{
  union
  {
    float scalars[vec_step(vector float)];
    vector float v;
  } buffer;
  buffer.scalars[0] = a;
  buffer.scalars[1] = b;
  buffer.scalars[2] = c;
  buffer.scalars[3] = d;
  return buffer.v;
}

void
GetVec (vector unsigned int i, int *i0, int *i1, int *i2, int *i3)
{
  union
  {
    unsigned int scalars[vec_step(vector float)];
    vector unsigned int v;
  } buffer;
  buffer.v = i;
  *i0 = buffer.scalars[0];
  *i1 = buffer.scalars[1];
  *i2 = buffer.scalars[2];
  *i3 = buffer.scalars[3];
}

vector unsigned char perm0 = (vector unsigned char)
                             ( 0, 1, 2, 3, 16, 17, 18, 19, 8, 9, 10, 11, 24, 25, 26, 27 );
vector unsigned char perm1 = (vector unsigned char)
                             (4, 5, 6, 7, 20, 21, 22, 23, 12, 13, 14, 15, 28, 29, 30, 31 );
vector unsigned char perm2 = (vector unsigned char)
                             ( 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 );
vector unsigned char perm3 = (vector unsigned char)
                             ( 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31 );
vector float zero = (vector float) (0.0, 0.0, 0.0, 0.0);

inline
vector float LoadUnaligned(float *target )
{
  vector float MSQ, LSQ, result;
  vector unsigned char mask;
  MSQ = vec_ld(0, target);          // most significant quadword
  LSQ = vec_ld(15, target);         // least significant quadword
  mask = vec_lvsl(0, target);       // create the permute mask
  result =  vec_perm(MSQ, LSQ, mask);  // align the data
  //  fprintf (stderr, "result = %vf\n", result);
  //   fprintf (stderr, "target = %f %f %f %f\n", target[0], target[1], target[2], target[3]);
  return result;
}



/// SSE3 add "horizontal add" instructions, which makes things
/// simpler and faster
// Use plain-old SSE instructions
#define _TRANSPOSE4(_v0, _v1, _v2, _v3)           \
do {                                              \
  vector float _t0 = vec_perm (_v0, _v1, perm0);  \
  vector float _t1 = vec_perm (_v0, _v1, perm1);  \
  vector float _t2 = vec_perm (_v2, _v3, perm0);  \
  vector float _t3 = vec_perm (_v2, _v3, perm1);  \
  _v0 = vec_perm (_t0, _t2, perm2);               \
  _v1 = vec_perm (_t1, _t3, perm2);               \
  _v2 = vec_perm (_t0, _t2, perm3);               \
  _v3 = vec_perm (_t1, _t3, perm3);               \
} while (0);

#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)                        \
do {                                                                \
  vector float r0 = vec_madd (M0, v, zero);                         \
  vector float r1 = vec_madd (M1, v, zero);              	    \
  vector float r2 = vec_madd (M2, v, zero);                         \
  vector float r3 = vec_madd (M3, v, zero);		            \
  _TRANSPOSE4 (r0, r1, r2, r3);                                     \
  r = vec_add (vec_add(r0, r1), vec_add(r2, r3));                   \
 } while (0);
#define _MM_DOT4_PS(A, B, p)                                        \
do {                                                                \
  vector float _t    = vec_madd (A, B, zero);                       \
  vector float _alo  = vec_mergel (_t, _t);                         \
  vector float _ahi  = vec_mergeh (_t, _t);                         \
  vector float _a    = vec_add (_alo, _ahi);                        \
  vector float _rlo  = vec_mergel (_a, _a);                         \
  vector float _rhi  = vec_mergeh (_a, _a);                         \
  vector float _r    = vec_add (_rlo, _rhi);                        \
  vector float _r2   = vec_splat (_r, 0);                           \
  vec_ste (_r2, 0, (p));                                            \
} while(0);

#define _4DOTS(_u0, _v0, _u1, _v1, _u2, _v2, _u3, _v3, result)      \
do {                                                                \
  vector float _w0   = vec_madd (_u0, _v0, zero);                   \
  vector float _w1   = vec_madd (_u1, _v1, zero);                   \
  vector float _w2   = vec_madd (_u2, _v2, zero);                   \
  vector float _w3   = vec_madd (_u3, _v3, zero);                   \
  _TRANSPOSE4 (_w0, _w1, _w2, _w3);                                 \
  result = vec_add (vec_add(_w0, _w1), vec_add(_w2, _w3));         \
} while(0);



/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_1d_s (UBspline_1d_s * restrict spline,
                    double x, float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;
  *val =
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
inline void
eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x,
                       float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;
  *val =
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv *
          (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
           coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
           coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
           coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}
/* Value, first derivative, and second derivative */
inline void
eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x,
                        float* restrict val, float* restrict grad,
                        float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  float* restrict coefs = spline->coefs;
  float tp[4];
  tp[0] = t*t*t;
  tp[1] = t*t;
  tp[2] = t;
  tp[3] = 1.0;
  *val =
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv *
          (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
           coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
           coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
           coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
          (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
           coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
           coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
           coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_2d_s (UBspline_2d_s * restrict spline,
                    double x, double y, float* restrict val)
{
}


/* Value and gradient */
inline void
eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline,
                       double x, double y,
                       float* restrict val, float* restrict grad)
{
}

/* Value, gradient, and laplacian */
inline void
eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline,
                        double x, double y, float* restrict val,
                        float* restrict grad, float* restrict lapl)
{
}

/* Value, gradient, and Hessian */
inline void
eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline,
                        double x, double y, float* restrict val,
                        float* restrict grad, float* restrict hess)
{
}



/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
inline void
eval_UBspline_3d_s (UBspline_3d_s * restrict spline,
                    double x, double y, double z,
                    float* restrict val)
{
}

/* Value and gradient */
inline void
eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline,
                       double x, double y, double z,
                       float* restrict val, float* restrict grad)
{
}



/* Value, gradient, and laplacian */
inline void
eval_UBspline_3d_s_vgl (UBspline_3d_s * restrict spline,
                        double x, double y, double z,
                        float* restrict val, float* restrict grad, float* restrict lapl)
{
}


/* Value, gradient, and Hessian */
inline void
eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline,
                        double x, double y, double z,
                        float* restrict val, float* restrict grad,
                        float* restrict hess)
{
  vec_dst (&A0, (12<<3) | (1<<8), 0);
  /// SSE mesh point determination
  vector float xyz       = MakeVec (x, y, z, 0.0);
  vector float x0y0z0    = MakeVec ( spline->x_grid.start,  spline->y_grid.start,
                                     spline->z_grid.start, 0.0);
  vector float delta_inv = MakeVec( spline->x_grid.delta_inv,
                                    spline->y_grid.delta_inv,
                                    spline->z_grid.delta_inv, 0.0 );
  xyz = vec_sub (xyz, x0y0z0);
  // ux = (x - x0)/delta_x and same for y and z
  vector float uxuyuz  = vec_madd (xyz, delta_inv, zero);
  //  fprintf (stderr, "uxuyuz = %vf\n", uxuyuz);
  // intpart = trunc (ux, uy, uz)
  vector float intpart  = vec_floor (uxuyuz);
  // fprintf (stderr, "intpart = %vf\n", intpart);
  vector unsigned int ixiyiz     = vec_ctu   (intpart, 0);
  // Store to memory for use in C expressions
  // xmm registers are stored to memory in reverse order
  int ix, iy, iz, dummy;
  //fprintf (stderr, "ixiyiz = %vld\n", ixiyiz);
  GetVec (ixiyiz, &ix, &iy, &iz, &dummy);
  // fprintf (stderr, "ix = %d  iy = %d  iz = %d\n", ix, iy, iz);
  int xs = spline->x_stride;
  int ys = spline->y_stride;
  // This macro is used to give the pointer to coefficient data.
  // i and j should be in the range [0,3].  Coefficients are read four
  // at a time, so no k value is needed.
#define P(i,j) ((float*)spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
  // Prefetch the data from main memory into cache so it's available
  // when we need to use it.
  int control_word;
  control_word = (2<<3) | (4<<8) | ((4*ys) << 16);
//   fprintf (stderr, "control word = %x\n", control_word);
//   fprintf (stderr, "ys = %d P(0,1)-P(0,0) = %d\n", ys,
//   P(0,1)-P(0,0));
  void *ptr = P(0,0);
  __dcbt (P(0,0), 0);
  __dcbt (P(0,1), 0);
  __dcbt (P(0,2), 0);
  __dcbt (P(0,3), 0);
  __dcbt (P(0,0), 12);
  __dcbt (P(0,1),12);
  __dcbt (P(0,2),12);
  __dcbt (P(0,3),12);
  __dcbt (P(1,0), 0);
  __dcbt (P(1,1), 0);
  __dcbt (P(1,2), 0);
  __dcbt (P(1,3), 0);
  __dcbt (P(1,0), 12);
  __dcbt (P(1,1),12);
  __dcbt (P(1,2),12);
  __dcbt (P(1,3),12);
  __dcbt (P(2,0), 0);
  __dcbt (P(2,1), 0);
  __dcbt (P(2,2), 0);
  __dcbt (P(2,3), 0);
  __dcbt (P(2,0), 12);
  __dcbt (P(2,1),12);
  __dcbt (P(2,2),12);
  __dcbt (P(2,3),12);
  __dcbt (P(3,0), 0);
  __dcbt (P(3,1), 0);
  __dcbt (P(3,2), 0);
  __dcbt (P(3,3), 0);
  __dcbt (P(3,0), 12);
  __dcbt (P(3,1),12);
  __dcbt (P(3,2),12);
  __dcbt (P(3,3),12);
//   vec_dstt (P(0,0), control_word, 0);
//   vec_dstt (P(1,0), control_word, 1);
//   vec_dstt (P(2,0), control_word, 2);
//   vec_dstt (P(3,0), control_word, 3);
//   // Now compute the vectors:
//   // tpx = [t_x^3 t_x^2 t_x 1]
//   // tpy = [t_y^3 t_y^2 t_y 1]
//   // tpz = [t_z^3 t_z^2 t_z 1]
  vector float txtytz = vec_sub (uxuyuz, intpart);
  vector float one    = (vector float) ( 1.0, 1.0, 1.0, 0.0);
  vector float t2     = vec_madd (txtytz, txtytz, zero);
  vector float t3     = vec_madd (t2, txtytz, zero);
//   vector float tpx    = t3;
//   vector float tpy    = t2;
//   vector float tpz    = txtytz;
//   vector float z2     = one;
//   _TRANSPOSE4(z2, tpz, tpy, tpx);
  vector float tpx    = t3;
  vector float tpy    = t2;
  vector float tpz    = txtytz;
  vector float z2     = one;
  _TRANSPOSE4(tpx, tpy, tpz, z2);
//   fprintf (stderr, "txtytz = %vf\n", txtytz);
//   fprintf (stderr, "tpxyz %vf   %vf   %vf\n", tpx, tpy, tpz);
//   fprintf (stderr, "ix,iy,iz = %d, %d, %d\n", ix, iy, iz);
  // a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
  // da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
  // A is 4x4 matrix given by the rows A0, A1, A2, A3
  vector float a, b, c, da, db, dc, d2a, d2b, d2c,
         cP[4], dcP[4], d2cP[4], bcP, dbcP, bdcP, d2bcP, dbdcP, bd2cP,
         tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  // x-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpx,   a);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpx,  da);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpx, d2a);
  // y-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpy,   b);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpy,  db);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpy, d2b);
  // z-dependent vectors
  _MM_MATVEC4_PS (  A0,   A1,   A2,   A3, tpz,   c);
  _MM_MATVEC4_PS ( dA0,  dA1,  dA2,  dA3, tpz,  dc);
  _MM_MATVEC4_PS (d2A0, d2A1, d2A2, d2A3, tpz, d2c);
//   fprintf (stderr, "a = %vf\n", a);
//   fprintf (stderr, "b = %vf\n", b);
//   fprintf (stderr, "c = %vf\n", c);
  // Compute cP, dcP, and d2cP products 1/4 at a time to maximize
  // register reuse and avoid rerereading from memory or cache.
  // 1st quarter
  tmp0 = LoadUnaligned (P(0,0));
  tmp1 = LoadUnaligned (P(0,1));
  tmp2 = LoadUnaligned (P(0,2));
  tmp3 = LoadUnaligned (P(0,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
  //  fprintf (stderr, "cP[0] = %vf\n", cP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[0]);
  // 2nd quarter
  tmp0 = LoadUnaligned (P(1,0));
  tmp1 = LoadUnaligned (P(1,1));
  tmp2 = LoadUnaligned (P(1,2));
  tmp3 = LoadUnaligned (P(1,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[1]);
  // 3rd quarter
  tmp0 = LoadUnaligned (P(2,0));
  tmp1 = LoadUnaligned (P(2,1));
  tmp2 = LoadUnaligned (P(2,2));
  tmp3 = LoadUnaligned (P(2,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[2]);
  // 4th quarter
  tmp0 = LoadUnaligned (P(3,0));
  tmp1 = LoadUnaligned (P(3,1));
  tmp2 = LoadUnaligned (P(3,2));
  tmp3 = LoadUnaligned (P(3,3));
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3, d2c, d2cP[3]);
  // Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);
  _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3], d2b, d2bcP);
  _MM_MATVEC4_PS (d2cP[0], d2cP[1], d2cP[2], d2cP[3],   b, bd2cP);
  _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],  db, dbdcP);
  vector float valgrad, hess4;
//   fprintf (stderr, "a = %vf\n", a);
//   fprintf (stderr, "bcP = %vf\n", bcP);
  _4DOTS (a, bcP, da, bcP, a, dbcP, a, bdcP, valgrad);
//   fprintf (stderr, "valgrad = %vf\n", valgrad);
  tmp0 = vec_splat (valgrad, 0);
  vec_ste (tmp0, 0, val);
  tmp0 = vec_splat (valgrad, 1);
  vec_ste (tmp0, 0, &(grad[0]));
  tmp0 = vec_splat (valgrad, 2);
  vec_ste (tmp0, 0, &(grad[1]));
  tmp0 = vec_splat (valgrad, 3);
  vec_ste (tmp0, 0, &(grad[2]));
  _4DOTS (d2a, bcP, a, d2bcP, a, bd2cP, da, dbcP, hess4);
  tmp0 = vec_splat (hess4, 0);
  vec_ste (tmp0, 0, &(hess[0]));
  tmp0 = vec_splat (hess4, 1);
  vec_ste (tmp0, 0, &(hess[4]));
  tmp0 = vec_splat (hess4, 2);
  vec_ste (tmp0, 0, &(hess[8]));
  tmp0 = vec_splat (hess4, 3);
  vec_ste (tmp0, 0, &(hess[1]));
  _4DOTS (da, bdcP, a, dbdcP, a, a, a, a, hess4);
  tmp0 = vec_splat (hess4, 0);
  vec_ste (tmp0, 0, &(hess[2]));
  tmp0 = vec_splat (hess4, 1);
  vec_ste (tmp0, 0, &(hess[5]));
  // Compute value
//   _MM_DOT4_PS (a, bcP, val);
//     // Compute gradient
//   _MM_DOT4_PS (da, bcP, &(grad[0]));
//   _MM_DOT4_PS (a, dbcP, &(grad[1]));
//   _MM_DOT4_PS (a, bdcP, &(grad[2]));
//   // Compute hessian
//   _MM_DOT4_PS (d2a, bcP, &(hess[0]));
//   _MM_DOT4_PS (a, d2bcP, &(hess[4]));
//   _MM_DOT4_PS (a, bd2cP, &(hess[8]));
//   _MM_DOT4_PS (da, dbcP, &(hess[1]));
//   _MM_DOT4_PS (da, bdcP, &(hess[2]));
//   _MM_DOT4_PS (a, dbdcP, &(hess[5]));
  // Multiply gradients and hessians by appropriate grid inverses
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;
  grad[0] *= dxInv;
  grad[1] *= dyInv;
  grad[2] *= dzInv;
  hess[0] *= dxInv*dxInv;
  hess[4] *= dyInv*dyInv;
  hess[8] *= dzInv*dzInv;
  hess[1] *= dxInv*dyInv;
  hess[2] *= dxInv*dzInv;
  hess[5] *= dyInv*dzInv;
  // Copy hessian elements into lower half of 3x3 matrix
  hess[3] = hess[1];
  hess[6] = hess[2];
  hess[7] = hess[5];
#undef P
  //fprintf (stderr, "%vf\n", xyz);
}

#undef _MM_MATVEC4_PS
#undef _MM_DOT4_PS

#endif
