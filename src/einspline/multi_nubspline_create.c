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

#include "multi_nubspline_create.h"
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#ifndef __USE_XOPEN2K
  #define __USE_XOPEN2K
#endif
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

int posix_memalign(void **memptr, size_t alignment, size_t size);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////       Helper functions for spline creation         ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void init_sse_data();


////////////////////////////////////////////////////////////
// Single-precision creation routines                     //
////////////////////////////////////////////////////////////
void
solve_NUB_deriv_interp_1d_s (NUBasis* restrict basis, 
			     float* restrict data, int datastride,
			     float* restrict    p, int pstride,
			     float abcdInitial[4], float abcdFinal[4]);
void
solve_NUB_periodic_interp_1d_s (NUBasis* restrict basis,
				float* restrict data, int datastride,
				float* restrict p, int pstride);

void
find_NUBcoefs_1d_s (NUBasis* restrict basis, BCtype_s bc,
		    float *data,  int dstride,
		    float *coefs, int cstride);


////////////////////////////////////////////////////////////
// Double-precision creation routines                     //
////////////////////////////////////////////////////////////
void
solve_NUB_deriv_interp_1d_d (NUBasis* restrict basis, 
			     double* restrict data, int datastride,
			     double* restrict    p, int pstride,
			     double abcdInitial[4], double abcdFinal[4]);

void
solve_NUB_periodic_interp_1d_d (NUBasis* restrict basis,
				double* restrict data, int datastride,
				double* restrict p, int pstride);

void
find_NUBcoefs_1d_d (NUBasis* restrict basis, BCtype_d bc,
		    double *data,  int dstride,
		    double *coefs, int cstride);



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////     Single-Precision, Real Creation Routines       ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
multi_NUBspline_1d_s*
create_multi_NUBspline_1d_s (NUgrid* x_grid, BCtype_s xBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_1d_s* restrict spline = malloc (sizeof(multi_NUBspline_1d_s));
  if (spline == NULL) 
    return spline;

  spline->spcode = MULTI_NU1D;
  spline->tcode  = SINGLE_REAL;
  
  // Next, create the basis
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->xBC = xBC; spline->x_grid = x_grid;
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx, Nx;
  if (xBC.lCode == PERIODIC)     Mx = x_grid->num_points - 1;
  else                           Mx = x_grid->num_points;
  Nx = x_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 4) 
    N += 4 - (N % 4);
#endif 

  spline->x_stride = N;
  spline->x_grid   = x_grid;
#ifndef HAVE_SSE
  spline->coefs = malloc (sizeof(float)*Nx*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, (sizeof(float)*Nx*N));
  init_sse_data();    
#endif

  return spline;
}

void
set_multi_NUBspline_1d_s (multi_NUBspline_1d_s *spline, int num,
			 float *data)
{
  float *coefs = spline->coefs + num;
  int xs = spline->x_stride;
  find_NUBcoefs_1d_s (spline->x_basis, spline->xBC, data, 1, 
		       coefs, xs);
}


multi_NUBspline_2d_s*
create_multi_NUBspline_2d_s (NUgrid* x_grid, NUgrid* y_grid,
			    BCtype_s xBC, BCtype_s yBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_2d_s* restrict spline = malloc (sizeof(multi_NUBspline_2d_s));
  spline->spcode = MULTI_NU2D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);

  int Mx, My, Nx, Ny;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 4) 
    N += 4 - (N % 4);
#endif

  spline->x_stride = Ny*N;
  spline->y_stride = N;
#ifndef HAVE_SSE
  spline->coefs = malloc ((size_t)sizeof(float)*Nx*Ny*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  sizeof(float)*Nx*Ny*N);
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_2d_s (multi_NUBspline_2d_s* spline, int num, float *data)
{
  int Mx, My, Nx, Ny;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;


  float *coefs = spline->coefs + num;
  int ys = spline->y_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy*ys;
    find_NUBcoefs_1d_s (spline->x_basis, spline->xBC, data+doffset, My,
			coefs+coffset, Ny*ys);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny*ys;
    intptr_t coffset = ix*Ny*ys;
    find_NUBcoefs_1d_s (spline->y_basis, spline->yBC, coefs+doffset, ys, 
			coefs+coffset, ys);
  }
}


multi_NUBspline_3d_s*
create_multi_NUBspline_3d_s (NUgrid* x_grid, NUgrid* y_grid, NUgrid* z_grid,
			    BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
			    int num_splines)
{
 // Create new spline
  multi_NUBspline_3d_s* restrict spline = malloc (sizeof(multi_NUBspline_3d_s));
  if (spline == NULL)
    return spline;
  spline->spcode = MULTI_NU3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->z_grid = z_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);
  spline->z_basis = create_NUBasis (z_grid, zBC.lCode==PERIODIC);

  int Mx, My, Mz, Nx, Ny, Nz;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;
  if (zBC.lCode == PERIODIC) Mz = z_grid->num_points - 1;
  else                       Mz = z_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;
  Nz = z_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 4) 
    N += 4 - (N % 4);
#endif

  spline->x_stride      = Ny*Nz*N;
  spline->y_stride      = Nz*N;
  spline->z_stride      = N;

#ifndef HAVE_SSE
  spline->coefs      = malloc (sizeof(float)*Nx*Ny*Nz*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  ((size_t)sizeof(float)*Nx*Ny*Nz*N));
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_3d_s (multi_NUBspline_3d_s* spline, int num, float *data)
{
  int Mx, My, Mz, Nx, Ny, Nz;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  if (spline->zBC.lCode == PERIODIC) Mz = spline->z_grid->num_points - 1;
  else                               Mz = spline->z_grid->num_points;

  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;
  Nz = spline->z_grid->num_points + 2;

  float *coefs = spline->coefs + num;

  int zs = spline->z_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      int doffset = iy*Mz+iz;
      int coffset = (iy*Nz+iz)*zs;
      find_NUBcoefs_1d_s (spline->x_basis, spline->xBC, data+doffset, My*Mz,
			  coefs+coffset, Ny*Nz*zs);
    }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      int doffset = (ix*Ny*Nz + iz)*zs;
      int coffset = (ix*Ny*Nz + iz)*zs;
      find_NUBcoefs_1d_s (spline->y_basis, spline->yBC, coefs+doffset, Nz*zs, 
			  coefs+coffset, Nz*zs);
    }

  // Now, solve in the Z-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      int doffset = ((ix*Ny+iy)*Nz)*zs;
      int coffset = ((ix*Ny+iy)*Nz)*zs;
      find_NUBcoefs_1d_s (spline->z_basis, spline->zBC, coefs+doffset, zs, 
			  coefs+coffset, zs);
    }
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////    Single-Precision, Complex Creation Routines     ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
multi_NUBspline_1d_c*
create_multi_NUBspline_1d_c (NUgrid* x_grid, BCtype_c xBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_1d_c* restrict spline = malloc (sizeof(multi_NUBspline_1d_c));
  if (spline == NULL)
    return spline;

  spline->spcode = MULTI_NU1D;
  spline->tcode  = SINGLE_COMPLEX;

  // Next, create the basis
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->xBC = xBC; 
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx, Nx;
  if (xBC.lCode == PERIODIC)     Mx = x_grid->num_points - 1;
  else                           Mx = x_grid->num_points;
  Nx = x_grid->num_points + 2;

  int N = num_splines;

#ifdef HAVE_SSE
  if (N % 2) 
    N += 2 - (N % 2);
#endif 

  spline->x_stride = N;
  spline->x_grid   = x_grid;

#ifndef HAVE_SSE
  spline->coefs = malloc (2*sizeof(float)*Nx*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 2*sizeof(float)*Nx*N);
  init_sse_data();    
#endif

  return spline;
}

void
set_multi_NUBspline_1d_c (multi_NUBspline_1d_c* spline, int num, 
			  complex_float *data)
{
  complex_float *coefs = spline->coefs + num;

  BCtype_s xBC_r, xBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;

  int xs = spline->x_stride;
  // Real part
  find_NUBcoefs_1d_s (spline->x_basis, xBC_r, 
		      (float*)data, 2, (float*)coefs, 2*xs);
  // Imaginarty part
  find_NUBcoefs_1d_s (spline->x_basis, xBC_i, 
		      ((float*)data)+1, 2, ((float*)coefs+1), 2*xs);
}



multi_NUBspline_2d_c*
create_multi_NUBspline_2d_c (NUgrid* x_grid, NUgrid* y_grid,
			    BCtype_c xBC, BCtype_c yBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_2d_c* restrict spline = malloc (sizeof(multi_NUBspline_2d_c));
  spline->spcode = MULTI_NU2D;
  spline->tcode  = SINGLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC;
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);

  // Setup internal variables
  int Mx, My, Nx, Ny;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 2)
    N++;
#endif

  spline->x_stride = Ny*N;
  spline->y_stride = N;

#ifndef HAVE_SSE
  spline->coefs = malloc (2*sizeof(float)*Nx*Ny*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  2*sizeof(float)*Nx*Ny*N);
#endif
  init_sse_data();

  return spline;
}


void
set_multi_NUBspline_2d_c (multi_NUBspline_2d_c* spline, int num, 
			  complex_float *data)
{
  // Setup internal variables
  int Mx, My, Nx, Ny;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;

  complex_float* coefs = spline->coefs + num;

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;
 
  int ys = spline->y_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    int doffset = (2*iy);
    int coffset = (2*iy)*ys;
    // Real part
    find_NUBcoefs_1d_s (spline->x_basis, xBC_r, ((float*)data)+doffset, 2*My,
			(float*)coefs+coffset, 2*Ny*ys);
    // Imag part
    find_NUBcoefs_1d_s (spline->x_basis, xBC_i, ((float*)data)+doffset+1, 2*My,
			((float*)coefs)+coffset+1, 2*Ny*ys);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    int doffset = (2*ix*Ny)*ys;
    int coffset = (2*ix*Ny)*ys;
    // Real part
    find_NUBcoefs_1d_s (spline->y_basis, yBC_r, ((float*)coefs)+doffset, 
			2*ys, ((float*)coefs)+coffset, 2*ys);
    // Imag part
    find_NUBcoefs_1d_s (spline->y_basis, yBC_i, ((float*)coefs)+doffset+1, 
			2*ys, ((float*)coefs)+coffset+1, 2*ys);
  }  
}

multi_NUBspline_3d_c*
create_multi_NUBspline_3d_c (NUgrid* x_grid, NUgrid* y_grid, NUgrid* z_grid,
		      BCtype_c xBC, BCtype_c yBC, BCtype_c zBC,
		      int num_splines)
{
  // Create new spline
  multi_NUBspline_3d_c* restrict spline = malloc (sizeof(multi_NUBspline_3d_c));
  spline->spcode = MULTI_NU3D;
  spline->tcode  = SINGLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->z_grid = z_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);
  spline->z_basis = create_NUBasis (z_grid, zBC.lCode==PERIODIC);

  int Mx, My, Mz, Nx, Ny, Nz;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;
  if (zBC.lCode == PERIODIC) Mz = z_grid->num_points - 1;
  else                       Mz = z_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;
  Nz = z_grid->num_points + 2;

  int N = spline->num_splines;
#ifdef HAVE_SSE
  if (N % 2)
    N++;
#endif

  spline->x_stride = Ny*Nz*N;
  spline->y_stride = Nz*N;
  spline->z_stride = N;

#ifndef HAVE_SSE
  spline->coefs = malloc ((size_t)2*sizeof(float)*Nx*Ny*Nz*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  (size_t)2*sizeof(float)*Nx*Ny*Nz*N);
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_3d_c (multi_NUBspline_3d_c* spline, int num, complex_float *data)
{
  int Mx, My, Mz, Nx, Ny, Nz;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  if (spline->zBC.lCode == PERIODIC) Mz = spline->z_grid->num_points - 1;
  else                               Mz = spline->z_grid->num_points;

  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;
  Nz = spline->z_grid->num_points + 2;

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;
  zBC_r.lCode = spline->zBC.lCode;  zBC_r.rCode = spline->zBC.rCode;
  zBC_r.lVal  = spline->zBC.lVal_r; zBC_r.rVal  = spline->zBC.rVal_r;
  zBC_i.lCode = spline->zBC.lCode;  zBC_i.rCode = spline->zBC.rCode;
  zBC_i.lVal  = spline->zBC.lVal_i; zBC_i.rVal  = spline->zBC.rVal_i;

  complex_float *coefs = spline->coefs + num;
  int zs = spline->z_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      int doffset = 2*(iy*Mz+iz);
      int coffset = 2*(iy*Nz+iz)*zs;
      // Real part
      find_NUBcoefs_1d_s (spline->x_basis, xBC_r, 
			  ((float*)data)+doffset, 2*My*Mz,
			  ((float*)coefs)+coffset, 2*Ny*Nz*zs);
      // Imag part
      find_NUBcoefs_1d_s (spline->x_basis, xBC_i, 
			  ((float*)data)+doffset+1, 2*My*Mz,
			  ((float*)coefs)+coffset+1, 2*Ny*Nz*zs);
    }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      int doffset = 2*(ix*Ny*Nz + iz)*zs;
      int coffset = 2*(ix*Ny*Nz + iz)*zs;
      // Real part
      find_NUBcoefs_1d_s (spline->y_basis, yBC_r, 
			  ((float*)coefs)+doffset, 2*Nz*zs, 
			  ((float*)coefs)+coffset, 2*Nz*zs);
      // Imag part
      find_NUBcoefs_1d_s (spline->y_basis, yBC_i, 
			  ((float*)coefs)+doffset+1, 2*Nz*zs, 
			  ((float*)coefs)+coffset+1, 2*Nz*zs);
    }

  // Now, solve in the Z-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      int doffset = 2*((ix*Ny+iy)*Nz)*zs;
      int coffset = 2*((ix*Ny+iy)*Nz)*zs;
      // Real part
      find_NUBcoefs_1d_s (spline->z_basis, zBC_r, 
			  ((float*)coefs)+doffset, 2*zs, 
			  ((float*)coefs)+coffset, 2*zs);
      // Imag part
      find_NUBcoefs_1d_s (spline->z_basis, zBC_i, 
			  ((float*)coefs)+doffset+1, 2*zs, 
			  ((float*)coefs)+coffset+1, 2*zs);
    }
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////     Double-Precision, Real Creation Routines       ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
multi_NUBspline_1d_d*
create_multi_NUBspline_1d_d (NUgrid* x_grid, BCtype_d xBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_1d_d* restrict spline = malloc (sizeof(multi_NUBspline_1d_d));
  if (spline == NULL)
    return spline;

  spline->spcode = MULTI_NU1D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->x_grid = x_grid;
  spline->num_splines = num_splines;

  // Next, create the basis
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);

  // Setup internal variables
  int Mx, Nx;
  if (xBC.lCode == PERIODIC)     Mx = x_grid->num_points - 1;
  else                           Mx = x_grid->num_points;
  Nx = x_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE2
  // We must pad to keep data aligned for SSE operations
  if (N & 1)
    N++;
#endif
  spline->x_stride = N;

#ifndef HAVE_SSE2
  spline->coefs = malloc (sizeof(double)*Nx*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, sizeof(double)*Nx*N);
  init_sse_data();
#endif
    
  return spline;
}

void
set_multi_NUBspline_1d_d (multi_NUBspline_1d_d* spline, int num, double *data)
{
  double *coefs = spline->coefs + num;
  int xs = spline->x_stride;
  find_NUBcoefs_1d_d (spline->x_basis, spline->xBC, data, 1, coefs, xs);
}

void
set_multi_NUBspline_1d_d_BC (multi_NUBspline_1d_d* spline, int num, double *data,
			     BCtype_d xBC)
{
  double *coefs = spline->coefs + num;
  int xs = spline->x_stride;
  find_NUBcoefs_1d_d (spline->x_basis, xBC, data, 1, coefs, xs);
}


multi_NUBspline_2d_d*
create_multi_NUBspline_2d_d (NUgrid* x_grid, NUgrid* y_grid,
			     BCtype_d xBC, BCtype_d yBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_2d_d* restrict spline = malloc (sizeof(multi_NUBspline_2d_d));
  spline->spcode = MULTI_NU2D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->num_splines = num_splines;
 
  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);

  int Mx, My, Nx, Ny;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE2
  // We must pad to keep data align for SSE operations
  if (num_splines & 1)
    N++;
#endif
  spline->x_stride = Ny*N;
  spline->y_stride = N;

#ifndef HAVE_SSE2
  spline->coefs = malloc (sizeof(double)*Nx*Ny*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, (sizeof(double)*Nx*Ny*N));
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_2d_d (multi_NUBspline_2d_d* spline, int num, double *data)
{
  int Mx, My, Nx, Ny;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;

  double *coefs = spline->coefs + num;
  int ys = spline->y_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    int doffset = iy;
    int coffset = iy*ys;
    find_NUBcoefs_1d_d (spline->x_basis, spline->xBC, data+doffset, My,
			coefs+coffset, Ny*ys);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    int doffset = ix*Ny*ys;
    int coffset = ix*Ny*ys;
    find_NUBcoefs_1d_d (spline->y_basis, spline->yBC, coefs+doffset, ys, 
		     coefs+coffset, ys);
  }
}


multi_NUBspline_3d_d*
create_multi_NUBspline_3d_d (NUgrid* x_grid, NUgrid* y_grid, NUgrid* z_grid,
			    BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
			    int num_splines)
{
  // Create new spline
  multi_NUBspline_3d_d* restrict spline = malloc (sizeof(multi_NUBspline_3d_d));
  if (spline == NULL)
    return spline;
  spline->spcode = MULTI_NU3D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->z_grid = z_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);
  spline->z_basis = create_NUBasis (z_grid, zBC.lCode==PERIODIC);

  int Mx, My, Mz, Nx, Ny, Nz;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;
  if (zBC.lCode == PERIODIC) Mz = z_grid->num_points - 1;
  else                       Mz = z_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;
  Nz = z_grid->num_points + 2;


  int N = num_splines;
#ifdef HAVE_SSE2
  // We must pad to keep data align for SSE operations
  if (N & 1)
    N++;
#endif
  
  spline->x_stride = Ny*Nz*N;
  spline->y_stride = Nz*N;
  spline->z_stride = N;
  
#ifndef HAVE_SSE2
  spline->coefs      = malloc ((size_t)sizeof(double)*Nx*Ny*Nz*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  ((size_t)sizeof(double)*Nx*Ny*Nz*N));
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_3d_d (multi_NUBspline_3d_d* spline, int num, double *data)
{
  int Mx, My, Mz, Nx, Ny, Nz;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  if (spline->zBC.lCode == PERIODIC) Mz = spline->z_grid->num_points - 1;
  else                               Mz = spline->z_grid->num_points;
  
  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;
  Nz = spline->z_grid->num_points + 2;
  
  double *coefs = spline->coefs + num;
  intptr_t zs = spline->z_stride;

  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      int doffset = iy*Mz+iz;
      int coffset = (iy*Nz+iz)*zs;
      find_NUBcoefs_1d_d (spline->x_basis, spline->xBC, data+doffset, My*Mz,
			  coefs+coffset, Ny*Nz*zs);
    }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      int doffset = (ix*Ny*Nz + iz)*zs;
      int coffset = (ix*Ny*Nz + iz)*zs;
      find_NUBcoefs_1d_d (spline->y_basis, spline->yBC, coefs+doffset, Nz*zs, 
			  coefs+coffset, Nz*zs);
    }

  // Now, solve in the Z-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      int doffset = (ix*Ny+iy)*Nz*zs;
      int coffset = (ix*Ny+iy)*Nz*zs;
      find_NUBcoefs_1d_d (spline->z_basis, spline->zBC, coefs+doffset, zs, 
			  coefs+coffset, zs);
    }
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////    Double-Precision, Complex Creation Routines     ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs


multi_NUBspline_1d_z*
create_multi_NUBspline_1d_z (NUgrid* x_grid, BCtype_z xBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_1d_z* restrict spline = malloc (sizeof(multi_NUBspline_1d_z));
  spline->spcode = MULTI_NU1D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC;
  spline->x_grid = x_grid;
  spline->num_splines = num_splines;

  // Next, create the basis
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);

  if (spline->x_basis->grid != x_grid) {
    fprintf (stderr, "Error in basis creation.\n");
    abort();
  }
  if (spline->x_basis == NULL) {
    fprintf (stderr, "Error creating basis in create_multi_NUBspline_1d_z.\n");
    abort();
  }
	     
  // Setup internal variables
  int Mx, Nx;
  if (xBC.lCode == PERIODIC)     Mx = x_grid->num_points - 1;
  else                           Mx = x_grid->num_points;
  Nx = x_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 2) 
    N ++;
#endif 

  spline->x_stride = N;
#ifndef HAVE_SSE2
  spline->coefs = malloc (2*sizeof(double)*Nx*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 2*sizeof(double)*Nx*N);
  init_sse_data();   
#endif

  return spline;
}

void
set_multi_NUBspline_1d_z (multi_NUBspline_1d_z* spline, int num, complex_double *data)
{
  complex_double *coefs = spline->coefs + num;

  BCtype_d xBC_r, xBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  int xs = spline->x_stride;
  // Real part
  find_NUBcoefs_1d_d (spline->x_basis, xBC_r, (double*)data, 2, 
		   ((double*)coefs),   2*xs);
  // Imaginary part
  find_NUBcoefs_1d_d (spline->x_basis, xBC_i, ((double*)data)+1, 2, 
		   ((double*)coefs)+1, 2*xs);
 
}

void
set_multi_NUBspline_1d_z_BC (multi_NUBspline_1d_z *spline, int num, 
			    complex_double *data, BCtype_z xBC)
{
  complex_double *coefs = spline->coefs + num;

  BCtype_d xBC_r, xBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  int xs = spline->x_stride;
  // Real part
  find_NUBcoefs_1d_d (spline->x_basis, xBC_r, (double*)data, 2, 
		      ((double*)coefs),   2*xs);
  // Imaginary part
  find_NUBcoefs_1d_d (spline->x_basis, xBC_i, ((double*)data)+1, 2, 
		      ((double*)coefs)+1, 2*xs);
}


multi_NUBspline_2d_z*
create_multi_NUBspline_2d_z (NUgrid* x_grid, NUgrid* y_grid,
			     BCtype_z xBC, BCtype_z yBC, int num_splines)
{
  // Create new spline
  multi_NUBspline_2d_z* restrict spline = malloc (sizeof(multi_NUBspline_2d_z));
  spline->spcode = MULTI_NU2D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC;
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);

  int Mx, My, Nx, Ny;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 4) 
    N += 4 - (N % 4);
#endif

  spline->x_stride = Ny*N;
  spline->y_stride = N;
  
#ifndef HAVE_SSE2
  spline->coefs = malloc (2*sizeof(double)*Nx*Ny*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 2*sizeof(double)*Nx*Ny*N);
  init_sse_data();
#endif

  return spline;
}


void
set_multi_NUBspline_2d_z (multi_NUBspline_2d_z* spline, int num,
			  complex_double *data)
{
  int Mx, My, Nx, Ny;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;

  complex_double *coefs = spline->coefs + num;
  int ys = spline->y_stride;

  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    int doffset = 2*iy;
    int coffset = 2*iy*ys;
    // Real part
    find_NUBcoefs_1d_d (spline->x_basis, xBC_r, 
			((double*)data+doffset), 2*My,
			(double*)coefs+coffset, 2*Ny*ys);
    // Imag part
    find_NUBcoefs_1d_d (spline->x_basis, xBC_i, ((double*)data)+doffset+1, 2*My,
			((double*)coefs)+coffset+1, 2*Ny*ys);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    int doffset = 2*ix*Ny*ys;
    int coffset = 2*ix*Ny*ys;
    // Real part
    find_NUBcoefs_1d_d (spline->y_basis, yBC_r, 
			((double*)coefs)+doffset, 2*ys, 
			(double*)coefs+coffset, 2*ys);
    // Imag part
    find_NUBcoefs_1d_d (spline->y_basis, yBC_i, 
			(double*)coefs+doffset+1, 2*ys, 
			((double*)coefs)+coffset+1, 2*ys);
  }
}



multi_NUBspline_3d_z*
create_multi_NUBspline_3d_z (NUgrid* x_grid, NUgrid* y_grid, NUgrid* z_grid,
			     BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
			     int num_splines)
{
  // Create new spline
  multi_NUBspline_3d_z* restrict spline = malloc (sizeof(multi_NUBspline_3d_z));
  spline->spcode = MULTI_NU3D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC;
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->z_grid = z_grid;
  spline->num_splines = num_splines;

  // Next, create the bases
  spline->x_basis = create_NUBasis (x_grid, xBC.lCode==PERIODIC);
  spline->y_basis = create_NUBasis (y_grid, yBC.lCode==PERIODIC);
  spline->z_basis = create_NUBasis (z_grid, zBC.lCode==PERIODIC);

  int Mx, My, Mz, Nx, Ny, Nz;
  if (xBC.lCode == PERIODIC) Mx = x_grid->num_points - 1;
  else                       Mx = x_grid->num_points;
  if (yBC.lCode == PERIODIC) My = y_grid->num_points - 1;
  else                       My = y_grid->num_points;
  if (zBC.lCode == PERIODIC) Mz = z_grid->num_points - 1;
  else                       Mz = z_grid->num_points;

  Nx = x_grid->num_points + 2;
  Ny = y_grid->num_points + 2;
  Nz = z_grid->num_points + 2;

  int N = num_splines;
#ifdef HAVE_SSE2
  if (N & 3)
    N += 4-(N & 3);
#endif

  spline->x_stride = Ny*Nz*N;
  spline->y_stride = Nz*N;
  spline->z_stride = N;

#ifndef HAVE_SSE2
  spline->coefs      = malloc ((size_t)2*sizeof(double)*Nx*Ny*Nz*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, (size_t)2*sizeof(double)*Nx*Ny*Nz*N);
  init_sse_data();
#endif

  return spline;
}

void
set_multi_NUBspline_3d_z (multi_NUBspline_3d_z* spline, int num, complex_double *data)
{
  int Mx, My, Mz, Nx, Ny, Nz;
  if (spline->xBC.lCode == PERIODIC) Mx = spline->x_grid->num_points - 1;
  else                               Mx = spline->x_grid->num_points;
  if (spline->yBC.lCode == PERIODIC) My = spline->y_grid->num_points - 1;
  else                               My = spline->y_grid->num_points;
  if (spline->zBC.lCode == PERIODIC) Mz = spline->z_grid->num_points - 1;
  else                               Mz = spline->z_grid->num_points;

  Nx = spline->x_grid->num_points + 2;
  Ny = spline->y_grid->num_points + 2;
  Nz = spline->z_grid->num_points + 2;

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;
  zBC_r.lCode = spline->zBC.lCode;  zBC_r.rCode = spline->zBC.rCode;
  zBC_r.lVal  = spline->zBC.lVal_r; zBC_r.rVal  = spline->zBC.rVal_r;
  zBC_i.lCode = spline->zBC.lCode;  zBC_i.rCode = spline->zBC.rCode;
  zBC_i.lVal  = spline->zBC.lVal_i; zBC_i.rVal  = spline->zBC.rVal_i;

  complex_double *coefs = spline->coefs + num;

  int N = spline->num_splines;
  int zs = spline->z_stride;
  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      int doffset = 2*(iy*Mz+iz);
      int coffset = 2*(iy*Nz+iz)*zs;
      // Real part
      find_NUBcoefs_1d_d (spline->x_basis, xBC_r, ((double*)data)+doffset, 2*My*Mz,
		       ((double*)coefs)+coffset, 2*Ny*Nz*zs);
      // Imag part
      find_NUBcoefs_1d_d (spline->x_basis, xBC_i, ((double*)data)+doffset+1, 2*My*Mz,
		       ((double*)coefs)+coffset+1, 2*Ny*Nz*zs);
    }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      int doffset = 2*(ix*Ny*Nz + iz)*zs;
      int coffset = 2*(ix*Ny*Nz + iz)*zs;
      // Real part
      find_NUBcoefs_1d_d (spline->y_basis, yBC_r, ((double*)coefs)+doffset, 2*Nz*zs, 
		       ((double*)coefs)+coffset, 2*Nz*zs);
      // Imag part
      find_NUBcoefs_1d_d (spline->y_basis, yBC_i, ((double*)coefs)+doffset+1, 2*Nz*zs, 
		       ((double*)coefs)+coffset+1, 2*Nz*zs);
    }

  // Now, solve in the Z-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      int doffset = 2*((ix*Ny+iy)*Nz)*zs;
      int coffset = 2*((ix*Ny+iy)*Nz)*zs;
      // Real part
      find_NUBcoefs_1d_d (spline->z_basis, zBC_r, ((double*)coefs)+doffset, 2*zs, 
		       ((double*)coefs)+coffset, 2*zs);
      // Imag part
      find_NUBcoefs_1d_d (spline->z_basis, zBC_i, ((double*)coefs)+doffset+1, 2*zs, 
		       ((double*)coefs)+coffset+1, 2*zs);
    }
}


void
destroy_multi_NUBspline (Bspline *spline)
{
  free (spline->coefs);
  free (spline);
}
