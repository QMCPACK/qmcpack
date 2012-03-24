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

#include "nubspline.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

double drand48();

void
PrintPassFail(bool pass)
{
  if (pass)
    // Print green "Passed"
    fprintf (stderr, "%c[32mPassed%c[0m\n", 0x1B, 0x1B);
  else
    // Print red "Failed"
    fprintf (stderr, "%c[31mFailed%c[0m\n", 0x1B, 0x1B);
}

void PrintTest (char *name, bool pass)
{
  int n = strlen (name);
  fprintf (stderr, "%s:", name);
  for (int i=n; i<57; i++)
    fprintf (stderr, " ");
  PrintPassFail (pass);
}


bool
TestCenterGrid()
{
  fprintf (stderr, "Testing CenterGrid:   ");
  bool passed = true;
  NUgrid* grid = create_center_grid (-5.0, 7.0, 6.0, 200);

  for (int i=0; i<10000; i++) {
    double x = -5.0+12.0*drand48();
    int lo = (*grid->reverse_map)(grid, x);
    assert (x >= grid->points[lo]);
    assert (x <= grid->points[lo+1]);
  }
  PrintPassFail (passed);
  return passed;
}


bool
TestGeneralGrid()
{
  fprintf (stderr, "Testing GeneralGrid:  ");
  bool passed = true;
  NUgrid* centgrid = create_center_grid (-5.0, 7.0, 6.0, 200);
  NUgrid* grid = create_general_grid (centgrid->points, 200);
  for (int i=0; i<10000; i++) {
    double x = -5.0+12.0*drand48();
    int lo = (*grid->reverse_map)(grid, x);
    passed = passed && (x >= grid->points[lo]);
    passed = passed && (x <= grid->points[lo+1]);
  }
  PrintPassFail (passed);
  return passed;
}

bool
close_float (float x, float y)
{
  float max = fmaxf (x, y);
  return (fabs(x-y)/max < 1.0e-5);
}

bool
TestNUB_1d_s()
{
  double start = -5.0;
  double end = 7.0;
  int N  = 200;
  NUgrid* grid = create_center_grid (start, end, 6.0, N);
  bool passed = true;
  float data[N];
  for (int i=0; i<N; i++) 
    data[i] = -1.0 + 2.0*drand48();
  BCtype_s bc;

  // Create spline with PBC
  fprintf (stderr, "Testing 1D single-precision periodic boundary conditions:\n");
  bc.lCode = PERIODIC; bc.rCode = PERIODIC;
  NUBspline_1d_s *periodic = create_NUBspline_1d_s (grid, bc, data);
  float sval, sgrad, slapl, eval, egrad, elapl;
  eval_NUBspline_1d_s_vgl (periodic, start, &sval, &sgrad, &slapl);
  eval_NUBspline_1d_s_vgl (periodic, end  , &eval, &egrad, &elapl);
  bool v_passed, grad_passed, lapl_passed;
  v_passed    = close_float (sval, eval);
  grad_passed = close_float (sgrad, egrad);
  lapl_passed = close_float (slapl, elapl);
  PrintTest ("Value", v_passed);
  PrintTest ("First derivative", grad_passed);
  PrintTest ("Second derivative", lapl_passed);
  passed = passed && v_passed && grad_passed && lapl_passed;

  double x = grid->points[26];
  float val;
  eval_NUBspline_1d_s (periodic, x, &val);
  bool interp_passed = close_float (val, data[26]);
  PrintTest ("Interpolation", interp_passed);
  passed = passed && interp_passed;

  // Create spline with fixed first derivative:
  bc.lCode = DERIV1; bc.lVal = 1.5;
  bc.rCode = DERIV1; bc.rVal = -0.3;
  NUBspline_1d_s *fixed_first = create_NUBspline_1d_s (grid, bc, data);
  fprintf (stderr, "Testing 1D single-precsion fixed first derivative boundary conditions:  \n");
  eval_NUBspline_1d_s_vg (fixed_first, start, &sval, &sgrad);
  eval_NUBspline_1d_s_vg (fixed_first,   end, &eval, &egrad);
  bool bc_passed = close_float (sgrad, 1.5) && close_float (egrad, -0.3);
  PrintTest ("Boundary conditions", bc_passed);
  x = grid->points[26];
  eval_NUBspline_1d_s (periodic, x, &val);
  interp_passed = close_float (val, data[26]);
  PrintTest ("Interpolation", interp_passed);
  passed = passed && interp_passed && bc_passed;

  // Create spline with fixed second derivative:
  bc.lCode = DERIV2; bc.lVal = 1.5;
  bc.rCode = DERIV2; bc.rVal = -0.3;
  NUBspline_1d_s *fixed_second = create_NUBspline_1d_s (grid, bc, data);
  fprintf (stderr, "Testing 1d_s fixed second derivative boundary conditions:  \n");
  eval_NUBspline_1d_s_vgl (fixed_second, start, &sval, &sgrad, &slapl);
  eval_NUBspline_1d_s_vgl (fixed_second,   end, &eval, &egrad, &elapl);
  bc_passed = close_float (slapl, 1.5) && close_float (elapl, -0.3);
  fprintf (stderr, "slapl = %1.8f  elapl = %1.8f\n", slapl, elapl);
  PrintTest ("Boundary conditions", bc_passed);
  x = grid->points[26];
  eval_NUBspline_1d_s (periodic, x, &val);
  interp_passed = close_float (val, data[26]);
  PrintTest ("Interpolation", interp_passed);
  passed = passed && interp_passed && bc_passed;

  return passed;
}

void
GridSpeedTest()
{
  NUgrid* centgrid = create_center_grid (-5.0, 7.0, 6.0, 2000);
  NUgrid* gengrid = create_general_grid (centgrid->points, 2000);
  int centsum=0, gensum=0;
  
  clock_t rstart, rend, cstart, cend, gstart, gend;
  
  rstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = -5.0 + 12.0*drand48();
  }
  rend = clock();

  cstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = -5.0 + 12.0*drand48();
    centsum += (*centgrid->reverse_map)(centgrid, x);
  }
  cend = clock();

  gstart = clock();
  for (int i=0; i<100000000; i++) {
    double x = -5.0 + 12.0*drand48();
    gensum += (*gengrid->reverse_map)(gengrid, x);
  }
  gend = clock();
  
  double cent_time = (double)(cend-cstart+rstart-rend)/(double)CLOCKS_PER_SEC;
  double gen_time  = (double)(gend-gstart+rstart-rend)/(double)CLOCKS_PER_SEC;
  fprintf (stderr, "%d %d\n", centsum, gensum);
  fprintf (stderr, "center_grid  time = %1.3f s.\n", cent_time);
  fprintf (stderr, "general_grid time = %1.3f s.\n", gen_time);
}

void
TestNUBasis()
{
  NUgrid* centgrid = create_center_grid (-5.0, 7.0, 10.0, 20);
  NUBasis* basis = create_NUBasis (centgrid, true);

  double bfuncs[4];
  for (double x=-5.0; x<=7.0; x+=0.001) {
    get_NUBasis_funcs_d (basis, x, bfuncs);
    fprintf (stderr, "%1.12f %1.12f %1.12f %1.12f %1.12f\n",
	     x, bfuncs[0], bfuncs[1], bfuncs[2], bfuncs[3]);
  }
}

void
TestNUBspline()
{
  NUgrid* centgrid = create_center_grid (-5.0, 7.0, 10.0, 20);
  NUBasis* basis = create_NUBasis (centgrid, true);
  float data[20];
  for (int i=0; i<20; i++) {
    double x = centgrid->points[i];
    double angle = (x+5.0)/12.0 * 2.0*M_PI;
    data[i] = sin(angle);
  }
  BCtype_s bc;
  //  bc.lCode = PERIODIC;  bc.rCode = PERIODIC;
  bc.lCode = DERIV1; bc.lVal = 2.0*M_PI/12.0;
  bc.rCode = DERIV1; bc.rVal = 2.0*M_PI/12.0;
  //bc.lCode = NATURAL;  bc.rCode = FLAT;
  NUBspline_1d_s *spline = create_NUBspline_1d_s (centgrid, bc, data);
  for (double x=-5.0; x<=7.0; x+=0.001) {
    float val, deriv;
    eval_NUBspline_1d_s_vg (spline, x, &val, &deriv);
    double angle = (x+5.0)/12.0 * 2.0*M_PI;
    fprintf (stderr, "%1.16e %1.16e %1.16e %1.16e\n", x, val, 
	     sin(angle), deriv);
  }
}


void
TestNUBspline_d()
{
  NUgrid* centgrid = create_center_grid (-5.0, 7.0, 10.0, 20);
  NUBasis* basis = create_NUBasis (centgrid, true);
  double data[20];
  for (int i=0; i<20; i++) {
    double x = centgrid->points[i];
    double angle = (x+5.0)/12.0 * 2.0*M_PI;
    data[i] = sin(angle);
  }
  BCtype_d bc;
  //  bc.lCode = PERIODIC;  bc.rCode = PERIODIC;
  bc.lCode = DERIV1; bc.lVal = 2.0*M_PI/12.0;
  bc.rCode = DERIV1; bc.rVal = 2.0*M_PI/12.0;
  //bc.lCode = NATURAL;  bc.rCode = FLAT;
  NUBspline_1d_d *spline = create_NUBspline_1d_d (centgrid, bc, data);
  for (double x=-5.0; x<=7.0; x+=0.001) {
    double val, deriv;
    eval_NUBspline_1d_d_vg (spline, x, &val, &deriv);
    double angle = (x+5.0)/12.0 * 2.0*M_PI;
    fprintf (stderr, "%1.16e %1.16e %1.16e %1.16e\n", x, val, 
	     sin(angle), deriv);
  }
}


void
TestNUB_2d_s()
{
  int Mx=30, My=35;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0, 7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0, 3.5, My);
  float data[Mx*My];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      data[ix*My+iy] = -1.0+2.0*drand48();
  
  BCtype_s xBC, yBC;
  xBC.lCode = PERIODIC;
  yBC.lCode = PERIODIC;
//   xBC.lCode = FLAT;  xBC.rCode = FLAT;
//   yBC.lCode = FLAT;  yBC.rCode = FLAT;

  NUBspline_2d_s *spline = create_NUBspline_2d_s (x_grid, y_grid, xBC, yBC, data);
  
  int xFine = 400;
  int yFine = 400;
  FILE *fout = fopen ("2d_s.dat", "w");
  double xi = x_grid->start;
  double xf = x_grid->end;// + x_grid->points[1] - x_grid->points[0];
  double yi = y_grid->start;
  double yf = y_grid->end;// + y_grid->points[1] - y_grid->points[0];
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      float val;
      eval_NUBspline_2d_s (spline, x, y, &val);
      fprintf (fout, "%1.16e ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
}


void
TestNUB_2d_c()
{
  int Mx=30, My=35;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0, 7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0, 3.5, My);
  complex_float data[Mx*My];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      data[ix*My+iy] = -1.0+2.0*drand48() + 1.0fi*(-1.0+2.0*drand48());
  
  BCtype_c xBC, yBC;
  xBC.lCode = PERIODIC;
  yBC.lCode = PERIODIC;
//   xBC.lCode = FLAT;  xBC.rCode = FLAT;
//   yBC.lCode = FLAT;  yBC.rCode = FLAT;

  NUBspline_2d_c *spline = create_NUBspline_2d_c (x_grid, y_grid, xBC, yBC, data);
  
  int xFine = 400;
  int yFine = 400;
  FILE *rout = fopen ("2d_r.dat", "w");
  FILE *iout = fopen ("2d_i.dat", "w");
  double xi = x_grid->start;
  double xf = x_grid->end;// + x_grid->points[1] - x_grid->points[0];
  double yi = y_grid->start;
  double yf = y_grid->end;// + y_grid->points[1] - y_grid->points[0];
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      complex_float val, grad[2], hess[4];
      eval_NUBspline_2d_c_vgh (spline, x, y, &val, grad, hess);
      fprintf (rout, "%1.16e ", crealf(val));
      fprintf (iout, "%1.16e ", cimagf(val));
    }
    fprintf (rout, "\n");
    fprintf (iout, "\n");
  }
  fclose (rout);
  fclose (iout);
}

void
TestNUB_3d_s()
{
  int Mx=20, My=27, Mz=23;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  3.5, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  2.8, Mz);
  float data[Mx*My*Mz];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48();
  
  BCtype_s xBC, yBC, zBC;
//   xBC.lCode = PERIODIC;
//   yBC.lCode = PERIODIC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_s *spline = create_NUBspline_3d_s (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  
  int xFine = 200, yFine = 200, zFine=200;
  FILE *fout = fopen ("3d_s.dat", "w");
  double xi = x_grid->start;  double xf = x_grid->end;
  double yi = y_grid->start;  double yf = y_grid->end;
  double zi = z_grid->start;  double zf = z_grid->end;
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      for (int iz=0; iz<zFine; iz++) {
	double z = zi + (double)iz/(double)(zFine)*(zf-zi);
	float val, grad[3], hess[9];
	eval_NUBspline_3d_s_vgh (spline, x, y, z, &val, grad, hess);
	fprintf (fout, "%1.16e ", val);
      }
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  fprintf (stderr, "spline->sp_code = %d\n", spline->sp_code);
  destroy_Bspline (spline);
}


void
TestNUB_3d_d()
{
  int Mx=20, My=27, Mz=23;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  3.5, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  2.8, Mz);
  double data[Mx*My*Mz];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48();
  
  BCtype_d xBC, yBC, zBC;
//   xBC.lCode = PERIODIC;
//   yBC.lCode = PERIODIC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_d *spline = create_NUBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  
  int xFine = 200, yFine = 200, zFine=200;
  FILE *fout = fopen ("3d_d.dat", "w");
  double xi = x_grid->start;  double xf = x_grid->end;
  double yi = y_grid->start;  double yf = y_grid->end;
  double zi = z_grid->start;  double zf = z_grid->end;
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      for (int iz=0; iz<zFine; iz++) {
	double z = zi + (double)iz/(double)(zFine)*(zf-zi);
	double val, grad[3], hess[9];
	eval_NUBspline_3d_d_vgh (spline, x, y, z, &val, grad, hess);
	fprintf (fout, "%1.16e ", val);
      }
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  fprintf (stderr, "spline->sp_code = %d\n", spline->sp_code);
  destroy_Bspline (spline);
}

void
TestNUB_3d_c()
{
  int Mx=20, My=27, Mz=23;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  3.5, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  2.8, Mz);
  complex_float data[Mx*My*Mz];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48() + 1.0if*(-1.0+2.0*drand48());
  
  BCtype_c xBC, yBC, zBC;
//   xBC.lCode = PERIODIC;
//   yBC.lCode = PERIODIC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_c *spline = create_NUBspline_3d_c (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  
  int xFine = 200, yFine = 200, zFine=200;
  FILE *rout = fopen ("3d_r.dat", "w");
  FILE *iout = fopen ("3d_i.dat", "w");
  double xi = x_grid->start;  double xf = x_grid->end;
  double yi = y_grid->start;  double yf = y_grid->end;
  double zi = z_grid->start;  double zf = z_grid->end;
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      for (int iz=0; iz<zFine; iz++) {
	double z = zi + (double)iz/(double)(zFine)*(zf-zi);
	complex_float val, grad[3], hess[9];
	eval_NUBspline_3d_c_vgh (spline, x, y, z, &val, grad, hess);
	fprintf (rout, "%1.16e ", crealf(val));
	fprintf (iout, "%1.16e ", cimagf(val));
      }
    }
    fprintf (rout, "\n");
    fprintf (iout, "\n");
  }
  fclose (rout);
  fclose (iout);
}


void
TestNUB_3d_z()
{
  int Mx=20, My=27, Mz=23;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  3.5, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  2.8, Mz);
  complex_double data[Mx*My*Mz];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48() + 1.0if*(-1.0+2.0*drand48());
  
  BCtype_z xBC, yBC, zBC;
//   xBC.lCode = PERIODIC;
//   yBC.lCode = PERIODIC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_z *spline = create_NUBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  
  int xFine = 200, yFine = 200, zFine=200;
  FILE *rout = fopen ("3d_r.dat", "w");
  FILE *iout = fopen ("3d_i.dat", "w");
  double xi = x_grid->start;  double xf = x_grid->end;
  double yi = y_grid->start;  double yf = y_grid->end;
  double zi = z_grid->start;  double zf = z_grid->end;
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      for (int iz=0; iz<zFine; iz++) {
	double z = zi + (double)iz/(double)(zFine)*(zf-zi);
	complex_double val, grad[3], hess[9];
	eval_NUBspline_3d_z_vgh (spline, x, y, z, &val, grad, hess);
	fprintf (rout, "%1.16e ", crealf(val));
	fprintf (iout, "%1.16e ", cimagf(val));
      }
    }
    fprintf (rout, "\n");
    fprintf (iout, "\n");
  }
  fclose (rout);
  fclose (iout);
}

void
SpeedNUB_3d_s()
{
  int Mx=200, My=200, Mz=200;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  1.0001, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  1.0001, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  1.0001, Mz);
  float *data;
  data = malloc (sizeof(float)*Mx*My*Mz);
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48();
  
  BCtype_s xBC, yBC, zBC;
//   xBC.lCode = PERIODIC;
//   yBC.lCode = PERIODIC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_s *spline = create_NUBspline_3d_s (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
 
  float val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid->start+ 0.9999*drand48()*(x_grid->end - x_grid->start);
    double y = y_grid->start+ 0.9999*drand48()*(y_grid->end - y_grid->start);
    double z = z_grid->start+ 0.9999*drand48()*(z_grid->end - z_grid->start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid->start+ 0.9999*drand48()*(x_grid->end - x_grid->start);
    double y = y_grid->start+ 0.9999*drand48()*(y_grid->end - y_grid->start);
    double z = z_grid->start+ 0.9999*drand48()*(z_grid->end - z_grid->start);
    eval_NUBspline_3d_s_vgh (spline, x, y, z, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
SpeedNUB_3d_z()
{
  int Mx=200, My=200, Mz=200;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0,  7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0,  3.5, My);
  NUgrid *z_grid = create_center_grid (-1.8, 2.0,  2.8, Mz);
  complex_double *data = malloc (sizeof(complex_double)*Mx*My*Mz);
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      for (int iz=0; iz<Mz; iz++)
	data[(ix*My+iy)*Mz+iz] = -1.0+2.0*drand48() + 1.0if*(-1.0+2.0*drand48());
  
  BCtype_z xBC, yBC, zBC;
  xBC.lCode = PERIODIC;  xBC.rCode = PERIODIC;
  yBC.lCode = PERIODIC;  yBC.rCode = PERIODIC;
  zBC.lCode = PERIODIC;  zBC.rCode = PERIODIC;

  NUBspline_3d_z *spline = create_NUBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, data);
  complex_double val, grad[3], hess[9];
  clock_t start, end, rstart, rend;
  rstart = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid->start+ 0.9999*drand48()*(x_grid->end - x_grid->start);
    double y = y_grid->start+ 0.9999*drand48()*(y_grid->end - y_grid->start);
    double z = z_grid->start+ 0.9999*drand48()*(z_grid->end - z_grid->start);
  }
  rend = clock();
  start = clock();
  for (int i=0; i<10000000; i++) {
    double x = x_grid->start+ 0.9999*drand48()*(x_grid->end - x_grid->start);
    double y = y_grid->start+ 0.9999*drand48()*(y_grid->end - y_grid->start);
    double z = z_grid->start+ 0.9999*drand48()*(z_grid->end - z_grid->start);
    eval_NUBspline_3d_z_vgh (spline, x, y, z, &val, grad, hess);
  }
  end = clock();
  fprintf (stderr, "10,000,000 evalations in %f seconds.\n", 
	   (double)(end-start-(rend-rstart))/(double)CLOCKS_PER_SEC);
}


void
TestNUB_2d_d()
{
  int Mx=30, My=35;
  NUgrid *x_grid = create_center_grid (-3.0, 4.0, 7.5, Mx);
  NUgrid *y_grid = create_center_grid (-1.0, 9.0, 3.5, My);
  double data[Mx*My];
  for (int ix=0; ix<Mx; ix++)
    for (int iy=0; iy<My; iy++)
      data[ix*My+iy] = -1.0+2.0*drand48();
  
  BCtype_d xBC, yBC;
  xBC.lCode = PERIODIC;
  yBC.lCode = PERIODIC;
//   xBC.lCode = FLAT;  xBC.rCode = FLAT;
//   yBC.lCode = FLAT;  yBC.rCode = FLAT;



  NUBspline_2d_d *spline = create_NUBspline_2d_d (x_grid, y_grid, xBC, yBC, data);
  
  int xFine = 400;
  int yFine = 400;
  FILE *fout = fopen ("2d_d.dat", "w");
  double xi = x_grid->start;
  double xf = x_grid->end;// + x_grid->points[1] - x_grid->points[0];
  double yi = y_grid->start;
  double yf = y_grid->end;// + y_grid->points[1] - y_grid->points[0];
  for (int ix=0; ix<xFine; ix++) {
    double x = xi+ (double)ix/(double)(xFine)*(xf-xi);
    for (int iy=0; iy<yFine; iy++) {
      double y = yi + (double)iy/(double)(yFine)*(yf-yi);
      double val;
      eval_NUBspline_2d_d (spline, x, y, &val);
      fprintf (fout, "%1.16e ", val);
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
}

int main()
{
  // TestCenterGrid();
  // TestGeneralGrid();
  // GridSpeedTest();
  // TestNUBasis();
  // TestNUBasis();
  TestNUBspline_d();
  // TestNUB_2d_s();
  //  TestNUB_2d_c();
  // TestNUB_3d_c();
  //  SpeedNUB_3d_s();
  // TestNUB_2d_d();
  // TestNUB_3d_d();
  // TestNUB_3d_z();
  //SpeedNUB_3d_z();
  //  bool passed = TestNUB_1d_s();
}

