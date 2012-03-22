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

#include "bspline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

double drand48();
void sincos (double phi, double *s, double *c);

typedef struct
{
  double kcut;
  double *Gvecs;
  double *coefs;
  int numG;
} periodic_func_d;

void
int_periodic_func (periodic_func_d *func, double kcut)
{
  func->kcut = kcut;
  func->numG = 0;
  int imax = (int) ceil (kcut/(2.0*M_PI));
  for (int ix=-imax; ix<=imax; ix++) {
    double kx = 2.0*M_PI * ix;
    for (int iy=-imax; iy<=imax; iy++) {
      double ky = 2.0*M_PI * iy;
      for (int iz=-imax; iz<=imax; iz++) {
	double kz = 2.0*M_PI * iz;
	if ((kx*kx + ky*ky + kz*kz) < (kcut*kcut))
	  func->numG++;
      }
    }
  }
  func->Gvecs = (double*) malloc (3*sizeof(double)*func->numG);
  func->coefs = (double*)  malloc (2*sizeof(double) *func->numG);

  int iG = 0;
  for (int ix=-imax; ix<=imax; ix++) {
    double kx = 2.0*M_PI * ix;
    for (int iy=-imax; iy<=imax; iy++) {
      double ky = 2.0*M_PI * iy;
      for (int iz=-imax; iz<=imax; iz++) {
	double kz = 2.0*M_PI * iz;
	if ((kx*kx + ky*ky + kz*kz) < (kcut*kcut)) {
	  func->Gvecs[3*iG+0] = kx;
	  func->Gvecs[3*iG+1] = ky;
	  func->Gvecs[3*iG+2] = kz;
	  func->coefs[2*iG+0] = 2.0*(drand48()-0.5);
	  func->coefs[2*iG+1] = 2.0*(drand48()-0.5);
	  iG++;
	}
      }
    }
  }
}

void
eval_periodic_func_d (periodic_func_d* restrict func,
		      double x, double y, double z,
		      double *restrict val, double *restrict grad,
		      double *restrict hess)
{
  *val = 0.0;
  for (int i=0; i<3; i++)    grad[i] = 0.0;
  for (int i=0; i<9; i++)    hess[i] = 0.0;

  for (int iG=0; iG<func->numG; iG++) {
    double kx = func->Gvecs[3*iG+0];
    double ky = func->Gvecs[3*iG+1];
    double kz = func->Gvecs[3*iG+2];
    double phase = x*kx + y*ky + z*kz;
    double re, im;
    sincos(phase, &im, &re);
    double c_re = func->coefs[2*iG+0];
    double c_im = func->coefs[2*iG+1];
    *val    += re*c_re - im*c_im;
    grad[0] += -kx*(re*c_im + im*c_re);
    grad[1] += -ky*(re*c_im + im*c_re);
    grad[2] += -kz*(re*c_im + im*c_re);
    hess[0] += -kx*kx*(re*c_re - im*c_im);
    hess[1] += -kx*ky*(re*c_re - im*c_im);
    hess[2] += -kx*kz*(re*c_re - im*c_im);
    hess[3] += -ky*kx*(re*c_re - im*c_im);
    hess[4] += -ky*ky*(re*c_re - im*c_im);
    hess[5] += -ky*kz*(re*c_re - im*c_im);
    hess[6] += -kz*kx*(re*c_re - im*c_im);
    hess[7] += -kz*ky*(re*c_re - im*c_im);
    hess[8] += -kz*kz*(re*c_re - im*c_im);
  }
}


void
test_bspline_3d_d()
{
  double kcut = 2.0*M_PI * 5.0;
  int Nspline = 100;
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num = Nspline;
  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num = Nspline;
  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num = Nspline;
  double dx = 1.0/(double)(Nspline);
  double dy = 1.0/(double)(Nspline);
  double dz = 1.0/(double)(Nspline);
  BCtype_d xBC, yBC, zBC;
  xBC.lCode = xBC.rCode = PERIODIC;
  yBC.lCode = yBC.rCode = PERIODIC;
  zBC.lCode = zBC.rCode = PERIODIC;

  double *data = malloc (sizeof(double)*Nspline*Nspline*Nspline);
  periodic_func_d func;
  int_periodic_func (&func, kcut);
  for (int ix=0; ix < x_grid.num; ix++) {
    double x = (double) ix * dx; 
    for (int iy=0; iy < y_grid.num; iy++) {
      double y = (double) iy * dy;
      for (int iz=0; iz < z_grid.num; iz++) {
	double z = (double) iz * dz;
	double val, grad[3], hess[9];
	eval_periodic_func_d (&func, x, y, z, &val, grad, hess);
	data[(ix*Nspline+iy)*Nspline+iz] = val;
      }
    }
  }
  
  UBspline_3d_d *spline =
    create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  
  int numTest = 10000;
  double valerror  = 0.0;
  double graderror = 0.0;
  double hesserror = 0.0;
  double valsum=0.0, gradsum=0.0, hesssum=0.0;
  for (int i=0; i<numTest; i++) {
    double x = drand48();
    double y = drand48();
    double z = drand48();
    double sval, sgrad[3], shess[9];
    double eval, egrad[3], ehess[9];
    
    eval_UBspline_3d_d_vgh (spline, x, y, z, &sval, sgrad, shess);
    eval_periodic_func_d   (&func,  x, y, z, &eval, egrad, ehess);
    valerror += (sval-eval)*(sval-eval);
    valsum += eval*eval;
    for (int i=0; i<3; i++) {
      graderror += (sgrad[i]-egrad[i])*(sgrad[i]-egrad[i]);
      gradsum   += egrad[i]*egrad[i];
    }
    for (int i=0; i<3; i++) {
      hesserror = (shess[i]-ehess[i])*(shess[i]-ehess[i]);
      hesssum += ehess[i]*ehess[i];
    }
    //    fprintf (stderr, "%10.8f %10.8f\n", eval, sval);
    //fprintf (stderr, "%14.8f %14.8f %14.8f     %14.8f %14.8f %14.8f\n",
    //	     egrad[0], egrad[1], egrad[2], sgrad[0], sgrad[1], sgrad[2]);
  }
  fprintf (stderr, "RMS val  error = %14.8f\n",
	   sqrt (valerror/valsum));
  fprintf (stderr, "RMS grad error = %14.8f\n",
	   sqrt (graderror/gradsum));
  fprintf (stderr, "RMS hess error = %14.8f\n",
	   sqrt (hesserror/hesssum));

}

main()
{
  test_bspline_3d_d();
}
