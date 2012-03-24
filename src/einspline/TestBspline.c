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

double drand48();

void
Test_1d_s()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 11;
  float data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  BCtype_s bc;
  bc.lCode = DERIV2; bc.lVal = 10.0;
  bc.rCode = DERIV2; bc.rVal = -10.0;
  
  FILE *fout = fopen ("1dSpline.dat", "w");
  UBspline_1d_s *spline = (UBspline_1d_s*) create_UBspline_1d_s (grid, bc, data);
  for (double x=1.0; x<=3.00001; x+=0.001) {
    float val, grad, lapl;
    eval_UBspline_1d_s_vgl (spline, x, &val, &grad, &lapl);
    fprintf (fout, "%1.5f %20.14f %20.14f %20.14f\n", x, val, grad, lapl);
  }
  fclose (fout);
}

void
Test_1d_d()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 1000;
  //  double data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  double data[10000];
  for (int i=0; i<10000; i++)
    data[i] = -2.0 + 4.0*drand48();
  BCtype_d bc;
  bc.lCode = DERIV1; bc.lVal = 10.0;
  bc.rCode = DERIV2; bc.rVal = -10.0;
  
  FILE *fout = fopen ("Spline_1d_d.dat", "w");
  UBspline_1d_d *spline = 
    (UBspline_1d_d*) create_UBspline_1d_d (grid, bc, data);
  for (double x=1.0; x<=3.00001; x+=0.001) {
    double val, grad, lapl;
    eval_UBspline_1d_d_vgl (spline, x, &val, &grad, &lapl);
    fprintf (fout, "%1.5f %20.14f %20.14f %20.14f\n", x, val, grad, lapl);
  }
  fclose (fout);
}

void
Test_1d_d_antiperiodic()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 10;
  //  double data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  double data[10];
  for (int i=0; i<10; i++)
    data[i] = -2.0 + 4.0*drand48();
  BCtype_d bc;
  bc.lCode = ANTIPERIODIC;
  
  FILE *fout = fopen ("Spline_1d_d_antiperiodic.dat", "w");
  UBspline_1d_d *spline = 
    (UBspline_1d_d*) create_UBspline_1d_d (grid, bc, data);
  for (double x=1.0; x<=5.00001; x+=0.001) {
    double val, grad, lapl;
    double xp = x;
    double sign = 1.0;
    while (xp >= grid.end) {
      xp -= (grid.end-grid.start);
      sign *= -1.0;
    }
    eval_UBspline_1d_d_vgl (spline, xp, &val, &grad, &lapl);
    fprintf (fout, "%1.5f %20.14f %20.14f %20.14f\n", x, sign*val, sign*grad, sign*lapl);
  }
  double val, grad, lapl;
  double x = grid.start + (grid.end-grid.start) * (double)1/(double)grid.num;
  eval_UBspline_1d_d_vgl (spline, x, &val, &grad, &lapl);
  fclose (fout);
}


void
Speed_1d_s()
{
  Ugrid grid;
  grid.start = 1.0;
  grid.end   = 3.0;
  grid.num = 11;
  float data[] = { 3.0, -4.0, 2.0, 1.0, -2.0, 0.0, 3.0, 2.0, 0.5, 1.0, 3.0 };
  BCtype_s bc;
  bc.lCode = DERIV2; bc.lVal = 10.0;
  bc.rCode = DERIV2; bc.rVal = -10.0;
  UBspline_1d_s *spline = (UBspline_1d_s*) create_UBspline_1d_s (grid, bc, data);

  float val, grad, lapl;
  clock_t start, end, rstart, rend;

  rstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = grid.start + 0.99999*drand48()*(grid.end-grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<100000000; i++) {
    double x = grid.start + 0.99999*drand48()*(grid.end-grid.start);
    eval_UBspline_1d_s_vgl (spline, x, &val, &grad, &lapl);
  }
  end = clock();
  fprintf (stderr, "100,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
Test_2d_s()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  
  float *data = malloc (x_grid.num * y_grid.num * sizeof(float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = -1.0 + 2.0*drand48();
  BCtype_s x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.lVal = 10.0;
  x_bc.rCode = PERIODIC; x_bc.rVal = -10.0;
  y_bc.lCode = PERIODIC; y_bc.lVal = 10.0;
  y_bc.rCode = PERIODIC; y_bc.rVal = -10.0;
  
  UBspline_2d_s *spline = (UBspline_2d_s*) create_UBspline_2d_s (x_grid, y_grid, x_bc, y_bc, data); 

  FILE *fout = fopen ("2dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      float val, grad[2], hess[4];
	eval_UBspline_2d_s_vgh (spline, x, y, &val, grad, hess);
      fprintf (fout, "%20.14f ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=5;
  int iy=7;
  float exval = data[ix*y_grid.num+iy];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  float spval, grad[2], hess[4];
  eval_UBspline_2d_s_vgh (spline, x, y, &spval, grad, hess);
  fprintf (stderr, "exval = %20.15f   spval = %20.15f\n", exval, spval);

}

void
Speed_2d_s()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 300;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 300;
  
  float *data = malloc (x_grid.num * y_grid.num * sizeof(float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = -1.0 + 2.0*drand48();
  BCtype_s x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC;
  
  UBspline_2d_s *spline = (UBspline_2d_s*) create_UBspline_2d_s (x_grid, y_grid, x_bc, y_bc, data); 
  float val, grad[2], hess[4];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    eval_UBspline_2d_s_vgh (spline, x, y, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}

void
Test_2d_c()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  
  complex_float *data = malloc (x_grid.num * y_grid.num * sizeof(complex_float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = 
	-1.0 + 2.0*drand48() + 1.0fI*(-1.0 + 2.0*drand48());
  BCtype_c x_bc, y_bc;
  x_bc.lCode = PERIODIC;  x_bc.rCode = PERIODIC;
  y_bc.lCode = PERIODIC;  y_bc.rCode = PERIODIC;
  
  UBspline_2d_c *spline = (UBspline_2d_c*) create_UBspline_2d_c (x_grid, y_grid, x_bc, y_bc, data); 

  FILE *fout = fopen ("2dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      complex_float val, grad[2], hess[4];
      eval_UBspline_2d_c_vgh (spline, x, y, &val, grad, hess);
      fprintf (fout, "%20.14f %20.15f ", crealf(val), cimagf(val));
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=5;
  int iy=7;
  complex_float exval = data[ix*y_grid.num+iy];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  complex_float spval, grad[2], hess[4];
  eval_UBspline_2d_c_vgh (spline, x, y, &spval, grad, hess);
  fprintf (stderr, "exval = (%20.15f + %20.15fi)   spval = (%20.15f + %20.15fi)\n", 
	   crealf(exval), cimagf(exval), creal(spval), cimagf(spval));

}

void
Speed_2d_c()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 300;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 300;
  
  complex_float *data = malloc (x_grid.num * y_grid.num * sizeof(complex_float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = 
	-1.0 + 2.0*drand48() + 1.0fI*(-1.0 + 2.0*drand48());
  BCtype_c x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC;
  
  UBspline_2d_c *spline = (UBspline_2d_c*) create_UBspline_2d_c (x_grid, y_grid, x_bc, y_bc, data); 
  complex_float val, grad[2], hess[4];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    eval_UBspline_2d_c_vgh (spline, x, y, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}

void
Test_2d_d()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  
  double *data = malloc (x_grid.num * y_grid.num * sizeof(double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = -1.0 + 2.0*drand48();
  BCtype_d x_bc, y_bc;
  x_bc.lCode = PERIODIC;  x_bc.rCode = PERIODIC;
  y_bc.lCode = PERIODIC;  y_bc.rCode = PERIODIC;
  
  UBspline_2d_d *spline = 
    create_UBspline_2d_d (x_grid, y_grid, x_bc, y_bc, data); 

  FILE *fout = fopen ("2dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      double val, grad[2], hess[4];
      eval_UBspline_2d_d_vgh (spline, x, y, &val, grad, hess);
      fprintf (fout, "%20.14f ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  
  int ix=5;
  int iy=7;
  double exval = data[ix*y_grid.num+iy];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  double spval, grad[2], hess[4];
  eval_UBspline_2d_d_vgh (spline, x, y, &spval, grad, hess);
  fprintf (stderr, "exval = %20.15f   spval = %20.15f\n", exval, spval);

}

void
Speed_2d_d()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 300;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 300;
  
  double *data = malloc (x_grid.num * y_grid.num * sizeof(double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = -1.0 + 2.0*drand48();
  BCtype_d x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC;
  
  UBspline_2d_d *spline = (UBspline_2d_d*) create_UBspline_2d_d (x_grid, y_grid, x_bc, y_bc, data); 
  double val, grad[2], hess[4];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<100000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    eval_UBspline_2d_d_vgh (spline, x, y, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
Test_2d_z()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  
  complex_double *data = malloc (x_grid.num * y_grid.num * sizeof(complex_double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = 
	-1.0 + 2.0*drand48() + 1.0I*(-1.0 + 2.0*drand48());
  BCtype_z x_bc, y_bc;
  x_bc.lCode = PERIODIC;  x_bc.rCode = PERIODIC;
  y_bc.lCode = PERIODIC;  y_bc.rCode = PERIODIC;
  
  UBspline_2d_z *spline = 
    create_UBspline_2d_z (x_grid, y_grid, x_bc, y_bc, data); 

  FILE *fout = fopen ("2dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      complex_double val, grad[2], hess[4];
      eval_UBspline_2d_z_vgh (spline, x, y, &val, grad, hess);
      fprintf (fout, "%20.14f %20.14f ", creal(val), cimag(val));
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  
  int ix=5;
  int iy=7;
  complex_double exval = data[ix*y_grid.num+iy];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  complex_double spval, grad[2], hess[4];
  eval_UBspline_2d_z_vgh (spline, x, y, &spval, grad, hess);
  fprintf (stderr, "exval = (%20.15f + %20.15fi)   spval = (%20.15f + %20.15fi)\n", 
	   creal(exval), cimag(exval), creal(spval), cimag(spval));

}

void
Speed_2d_z()
{
  Ugrid x_grid, y_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 300;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 300;
  
  complex_double *data = malloc (x_grid.num * y_grid.num * sizeof(complex_double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      *(data + ix*y_grid.num + iy) = 
	-1.0 + 2.0*drand48() + 1.0I*(-1.0 + 2.0*drand48());
  BCtype_z x_bc, y_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC;
  
  UBspline_2d_z *spline = (UBspline_2d_z*) create_UBspline_2d_z (x_grid, y_grid, x_bc, y_bc, data); 
  complex_double val, grad[2], hess[4];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<100000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    eval_UBspline_2d_z_vgh (spline, x, y, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "100,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}



void
Test_3d_s()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0001;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0001;  y_grid.num = 30;
  z_grid.start = 1.0;  z_grid.end   = 3.0001;  z_grid.num = 30;
  
  float *data = malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = -1.0 + 2.0*drand48();
  BCtype_s x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_s *spline = (UBspline_3d_s*) create_UBspline_3d_s 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  double z = 1.92341;
  FILE *fout = fopen ("3dspline.dat", "w");
  for (double x=x_grid.start; x<x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<y_grid.end; y+=0.005) {
      float val, grad[3], hess[9], lapl;
      eval_UBspline_3d_s_vgh (spline, x, y, z, &val, grad, hess);
      fprintf (fout, "%20.14f ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=9;  int iy=19; int iz = 24;
  float exval = data[(ix*y_grid.num+iy)*z_grid.num+iz];
  double x = x_grid.start + (double)ix * spline->x_grid.delta + 0.000001;
  double y = y_grid.start + (double)iy * spline->y_grid.delta + 0.000001;
  z =        z_grid.start + (double)iz * spline->z_grid.delta + 0.000001;
  float spval, grad[3], hess[9], lapl;
  eval_UBspline_3d_s_vgh (spline, x, y, z, &spval, grad, hess);
  fprintf (stderr, "exval = %20.15f   spval = %20.15f\n", exval, spval);

}


void
Speed_3d_s()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 200;
  y_grid.start = 1.0;  y_grid.end   = 5.0;  y_grid.num = 200;
  z_grid.start = 1.0;  z_grid.end   = 7.0;  z_grid.num = 200;
  
  float *data = malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = -1.0 + 2.0*drand48();
  BCtype_s x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_s *spline = (UBspline_3d_s*) create_UBspline_3d_s 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  float val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
    eval_UBspline_3d_s_vgh (spline, x, y, z, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
Test_3d_d()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  z_grid.start = 1.0;  z_grid.end   = 3.0;  z_grid.num = 30;
  
  double *data = malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = -1.0 + 2.0*drand48();
  BCtype_d x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_d *spline = (UBspline_3d_d*) create_UBspline_3d_d 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  double z = 1.92341;
  FILE *fout = fopen ("3dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      double val, grad[3], hess[9];
      eval_UBspline_3d_d_vgh (spline, x, y, z, &val, grad, hess);
      fprintf (fout, "%23.17f ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=9;  int iy=19; int iz = 24;
  double exval = data[(ix*y_grid.num+iy)*z_grid.num+iz];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  z =        z_grid.start + (double)iz * spline->z_grid.delta;
  double spval, grad[3], hess[9];
  eval_UBspline_3d_d_vgh (spline, x, y, z, &spval, grad, hess);
  fprintf (stderr, "exval = %23.17f   spval = %23.17f\n", exval, spval);

}


void
Speed_3d_d()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 200;
  y_grid.start = 1.0;  y_grid.end   = 5.0;  y_grid.num = 200;
  z_grid.start = 1.0;  z_grid.end   = 7.0;  z_grid.num = 200;
  
  double *data = malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = -1.0 + 2.0*drand48();
  BCtype_d x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_d *spline = (UBspline_3d_d*) create_UBspline_3d_d 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  double val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
    eval_UBspline_3d_d_vgh (spline, x, y, z, &val, grad, hess);
    // eval_UBspline_3d_d (spline, x, y, z, &val);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
Test_3d_c()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.0;  y_grid.num = 30;
  z_grid.start = 1.0;  z_grid.end   = 3.0;  z_grid.num = 30;
  
  complex_float *data = 
    malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(complex_float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = 
	  (-1.0 + 2.0*drand48()) + (-1.0 + 2.0*drand48())*1.0fI;
  BCtype_c x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_c *spline = create_UBspline_3d_c 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

   double z = 1.92341; 
  FILE *fout = fopen ("3dspline.dat", "w");
  for (double x=x_grid.start; x<0.99999*x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<y_grid.end; y+=0.005) {
      complex_float val, grad[3], hess[9];
      eval_UBspline_3d_c_vgh (spline, x, y, z, &val, grad, hess);
      fprintf (fout, "%23.17f %23.17f ", crealf(val), cimagf(val));
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=9;  int iy=18; int iz = 24;
  complex_float exval = data[(ix*y_grid.num+iy)*z_grid.num+iz];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  z =        z_grid.start + (double)iz * spline->z_grid.delta;
  complex_float spval, grad[3], hess[9];
  eval_UBspline_3d_c_vgh (spline, x, y, z, &spval, grad, hess);
  fprintf (stderr, "exval = (%23.17f + %23.17fi)\nspval = (%23.17f + %23.17fi)\n", 
	   crealf(exval), cimagf(exval), crealf(spval), cimagf(spval));

}


void
Speed_3d_c()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 200;
  y_grid.start = 1.0;  y_grid.end   = 5.0;  y_grid.num = 200;
  z_grid.start = 1.0;  z_grid.end   = 7.0;  z_grid.num = 200;
  
  complex_float *data = malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(complex_float));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = 
	  (-1.0 + 2.0*drand48()) + (-1.0 + 2.0*drand48())*1.0fI;
  BCtype_c x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_c *spline = (UBspline_3d_c*) create_UBspline_3d_c 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  complex_float val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
    eval_UBspline_3d_c_vgh (spline, x, y, z, &val, grad, hess);
    //eval_UBspline_3d_c     (spline, x, y, z, &val);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
Test_3d_z()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.4;  x_grid.num = 30;
  y_grid.start = 1.0;  y_grid.end   = 3.7;  y_grid.num = 30;
  z_grid.start = 1.0;  z_grid.end   = 3.9;  z_grid.num = 30;
  
  complex_double *data = 
    malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(complex_double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = 
	  (-1.0 + 2.0*drand48()) + (-1.0 + 2.0*drand48())*1.0fI;
  BCtype_z x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_z *spline = create_UBspline_3d_z 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  double z = 1.92341;
  FILE *fout = fopen ("3dspline.dat", "w");
  for (double x=x_grid.start; x<=x_grid.end; x+=0.005) {
    for (double y=y_grid.start; y<=y_grid.end; y+=0.005) {
      complex_double val, grad[3], hess[9];
      eval_UBspline_3d_z_vgh (spline, x, y, z, &val, grad, hess);
      fprintf (fout, "%23.19f %23.19f ", crealf(hess[4]), cimagf(hess[4]));
    }
    fprintf (fout, "\n");
  }
  fclose (fout);

  int ix=9;  int iy=19; int iz = 25;
  complex_double exval = data[(ix*y_grid.num+iy)*z_grid.num+iz];
  double x = x_grid.start + (double)ix * spline->x_grid.delta;
  double y = y_grid.start + (double)iy * spline->y_grid.delta;
  z =        z_grid.start + (double)iz * spline->z_grid.delta;
  complex_double spval, grad[3], hess[9];
  eval_UBspline_3d_z_vgh (spline, x, y, z, &spval, grad, hess);
  fprintf (stderr, "exval = (%23.19f + %23.19fi)\nspval = (%23.17f + %23.17fi)\n", 
	   crealf(exval), cimagf(exval), crealf(spval), cimagf(spval));

}


void
Speed_3d_z()
{
  Ugrid x_grid, y_grid, z_grid;
  x_grid.start = 1.0;  x_grid.end   = 3.0;  x_grid.num = 200;
  y_grid.start = 1.0;  y_grid.end   = 5.0;  y_grid.num = 200;
  z_grid.start = 1.0;  z_grid.end   = 7.0;  z_grid.num = 200;
  
  complex_double *data = 
    malloc (x_grid.num * y_grid.num * z_grid.num * sizeof(complex_double));
  for (int ix=0; ix<x_grid.num; ix++)
    for (int iy=0; iy<y_grid.num; iy++)
      for (int iz=0; iz<z_grid.num; iz++)
	*(data + ((ix*y_grid.num) + iy)*z_grid.num + iz) = 
	  (-1.0 + 2.0*drand48()) + (-1.0 + 2.0*drand48())*1.0fI;
  BCtype_z x_bc, y_bc, z_bc;
  x_bc.lCode = PERIODIC; x_bc.rCode = PERIODIC; 
  y_bc.lCode = PERIODIC; y_bc.rCode = PERIODIC; 
  z_bc.lCode = PERIODIC; z_bc.rCode = PERIODIC; 
  
  UBspline_3d_z *spline = (UBspline_3d_z*) create_UBspline_3d_z 
    (x_grid, y_grid, z_grid, x_bc, y_bc, z_bc, data); 

  complex_double val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid.start+ 0.9999*drand48()*(x_grid.end - x_grid.start);
    double y = y_grid.start+ 0.9999*drand48()*(y_grid.end - y_grid.start);
    double z = z_grid.start+ 0.9999*drand48()*(z_grid.end - z_grid.start);
    eval_UBspline_3d_z_vgh (spline, x, y, z, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}

#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
     extern "C"
#  endif
   int F77_DUMMY_MAIN() { return 1; }
#endif

int main()
{
  Test_1d_s();
  Test_1d_d();
  Test_1d_d_antiperiodic();
  // Speed_1d_s();
  Test_2d_s();
  // Speed_2d_s();
  Test_2d_c();
  // Speed_2d_c();
  Test_2d_d();
  // Speed_2d_d();
   Test_2d_z();
  // Speed_2d_z();
  Test_3d_s();
  // Speed_3d_s();
  Test_3d_d();
  // Speed_3d_d();
  Test_3d_c();
  // Speed_3d_c();
  Test_3d_z();
  Speed_3d_z();
}
