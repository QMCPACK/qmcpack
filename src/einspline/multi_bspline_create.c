/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include "multi_bspline_create.h"
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#ifndef __USE_XOPEN2K
  #define __USE_XOPEN2K
#endif
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

int posix_memalign(void **memptr, size_t alignment, size_t size);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////       Helper functions for spline creation         ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void init_sse_data();

void
find_coefs_1d_d (Ugrid grid, BCtype_d bc, 
		 double *data,  intptr_t dstride,
		 double *coefs, intptr_t cstride);

void 
solve_deriv_interp_1d_s (float bands[], float coefs[],
			 int M, int cstride);

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_periodic_interp_1d_s (float bands[], float coefs[],
			    int M, int cstride);

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_antiperiodic_interp_1d_s (float bands[], float coefs[],
				int M, int cstride);

void
find_coefs_1d_s (Ugrid grid, BCtype_s bc, 
		 float *data,  intptr_t dstride,
		 float *coefs, intptr_t cstride);

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
multi_UBspline_1d_s*
create_multi_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, int num_splines)
{
  // Create new spline
  multi_UBspline_1d_s* restrict spline = malloc (sizeof(multi_UBspline_1d_s));
  if (!spline) {
    fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_1d_s.\n");
    abort();
  }
  spline->spcode = MULTI_U1D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; spline->x_grid = x_grid;
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx = x_grid.num;
  int Nx;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    Nx = Mx+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    Nx = Mx+2;
  }

  int N = num_splines;
#ifdef HAVE_SSE
  if (N % 4) 
    N += 4 - (N % 4);
#endif 

  spline->x_stride = N;
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
#ifndef HAVE_POSIX_MEMALIGN
  spline->coefs = malloc (sizeof(float)*Nx*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, (sizeof(float)*Nx*N));
#endif
  spline->coefs_size=(size_t)Nx*(size_t)N;
#ifdef HAVE_SSE
  init_sse_data();    
#endif
  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficient in create_multi_UBspline_1d_s.\n");
    abort();
  }


  return spline;
}







multi_UBspline_3d_s*
create_multi_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
			    BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
			    int num_splines)
{
  // Create new spline
  multi_UBspline_3d_s* restrict spline = malloc (sizeof(multi_UBspline_3d_s));
  if (!spline) {
    fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_3d_s.\n");
    abort();
  }
  spline->spcode = MULTI_U3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->num_splines = num_splines;
  // Setup internal variables
  int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                           
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                           
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                           
    Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  int N = num_splines;
#if defined(HAVE_SSE)
  if (N % 4) 
    N += 4 - (N % 4);
  // fprintf(stdout, " The coefs has been 16-byte aligned.\n");
#endif

  spline->x_stride      = Ny*Nz*N;
  spline->y_stride      = Nz*N;
  spline->z_stride      = N;

#ifndef HAVE_POSIX_MEMALIGN
  spline->coefs      = malloc (sizeof(float)*Nx*Ny*Nz*N);
#else
  posix_memalign ((void**)&spline->coefs, 64, 
		  ((size_t)sizeof(float)*Nx*Ny*Nz*N));
#endif
  spline->coefs_size=(size_t)Nx*(size_t)Ny*(size_t)Nz*(size_t)N;
#ifdef HAVE_SSE
  init_sse_data();
#endif
  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_3d_s.\n");
    abort();
  }

  return spline;
}

void
set_multi_UBspline_3d_s (multi_UBspline_3d_s* spline, int num, float *data)
{
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                                   
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                                   
    Ny = My+2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                                   
    Nz = Mz+2;

  float *coefs = spline->coefs + num;

  intptr_t zs = spline->z_stride;
  // First, solve in the X-direction 
#pragma omp parallel for
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = (iy*Nz+iz)*zs;
      find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, (intptr_t)(My*Mz),
		       coefs+coffset, (intptr_t)(Ny*Nz)*zs);
    }
  
  // Now, solve in the Y-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      intptr_t doffset = (ix*Ny*Nz + iz)*zs;
      intptr_t coffset = (ix*Ny*Nz + iz)*zs;
      find_coefs_1d_s (spline->y_grid, spline->yBC, coefs+doffset, (intptr_t)Nz*zs, 
		       coefs+coffset, (intptr_t)Nz*zs);
    }

  // Now, solve in the Z-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      intptr_t doffset = ((ix*Ny+iy)*Nz)*zs;
      intptr_t coffset = ((ix*Ny+iy)*Nz)*zs;
      find_coefs_1d_s (spline->z_grid, spline->zBC, coefs+doffset, 
		       zs, coefs+coffset, zs);
    }
}


void
set_multi_UBspline_3d_s_d(multi_UBspline_3d_s* spline, int num, double *data)
{
  
  BCtype_d xBC, yBC, zBC;
  xBC.lCode=spline->xBC.lCode; xBC.rCode=spline->xBC.rCode;
  yBC.lCode=spline->yBC.lCode; yBC.rCode=spline->yBC.rCode;
  zBC.lCode=spline->zBC.lCode; zBC.rCode=spline->zBC.rCode;
  xBC.lVal=spline->xBC.lVal; xBC.rVal=spline->xBC.rVal;
  yBC.lVal=spline->yBC.lVal; yBC.rVal=spline->yBC.rVal;
  zBC.lVal=spline->zBC.lVal; zBC.rVal=spline->zBC.rVal;

  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)
    Nx = Mx+3;
  else
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)
    Ny = My+3;
  else
    Ny = My+2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)
    Nz = Mz+3;
  else
    Nz = Mz+2;

  double *spline_tmp = malloc(sizeof(double)*Nx*Ny*Nz);

  // First, solve in the X-direction 
#pragma omp parallel for
  for (int iy=0; iy<My; iy++)
    for (int iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_d (spline->x_grid, xBC, data+doffset, My*Mz, spline_tmp+coffset, Ny*Nz);
    }

  // Now, solve in the Y-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++)
    for (int iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_d (spline->y_grid, yBC, spline_tmp+doffset, Nz, spline_tmp+coffset, Nz);
    }

  // Now, solve in the Z-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_d (spline->z_grid, zBC, spline_tmp+doffset, 1, spline_tmp+coffset, 1);
    }

  {
//    const double* restrict i_ptr=spline_tmp;
#pragma omp parallel for
    for(int ix=0; ix<Nx; ++ix)
    {
      const double* restrict i_ptr=spline_tmp+ix*Ny*Nz;
      for(int iy=0; iy<Ny; ++iy)
        for(int iz=0; iz<Nz; ++iz)
          spline->coefs[ix*spline->x_stride +
                        iy*spline->y_stride +
                        iz*spline->z_stride + num] = (float)(*i_ptr++);
    }
  }

 free (spline_tmp);
}


/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////    Single-Precision, Complex Creation Routines     ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////     Double-Precision, Real Creation Routines       ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_deriv_interp_1d_d (double bands[], double coefs[],
			 int M, int cstride);

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void 
solve_periodic_interp_1d_d (double bands[], double coefs[],
			    int M, intptr_t cstride);

void
find_coefs_1d_d (Ugrid grid, BCtype_d bc, 
		 double *data,  intptr_t dstride,
		 double *coefs, intptr_t cstride);

multi_UBspline_1d_d*
create_multi_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, int num_splines)
{
  // Create new spline
  multi_UBspline_1d_d* restrict spline = malloc (sizeof(multi_UBspline_1d_d));
  if (!spline) {
    fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_1d_d.\n");
    abort();
  }
  spline->spcode = MULTI_U1D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx = x_grid.num;
  int Nx;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    Nx = Mx+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    Nx = Mx+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  int N = num_splines;
#ifdef HAVE_SSE2
  // We must pad to keep data aligned for SSE operations
  if (N & 1)
    N++;
#endif
  spline->x_stride = N;

#ifndef HAVE_POSIX_MEMALIGN
  spline->coefs = malloc (sizeof(double)*Nx*N);

#else
  posix_memalign ((void**)&spline->coefs, 64, sizeof(double)*Nx*N);
#endif
  spline->coefs_size=(size_t)Nx*(size_t)N;
#ifdef HAVE_SSE2
  init_sse_data();
#endif
  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_1d_d.\n");
    abort();
  }

  return spline;
}

void
set_multi_UBspline_1d_d (multi_UBspline_1d_d* spline, int num, double *data)
{
  double *coefs = spline->coefs + num;
  int xs = spline->x_stride;
  find_coefs_1d_d (spline->x_grid, spline->xBC, data, 1, coefs, xs);
}







multi_UBspline_3d_d*
create_multi_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
			    BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
			    int num_splines)
{
  // Create new spline
  multi_UBspline_3d_d* restrict spline;
#ifdef HAVE_POSIX_MEMALIGN
  posix_memalign ((void**)&spline, 64, (size_t)sizeof(multi_UBspline_3d_d));
#else
  spline = malloc (sizeof(multi_UBspline_3d_d));
#endif
  if (!spline) {
    fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_3d_d.\n");
    abort();
  }
  spline->spcode = MULTI_U3D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->num_splines = num_splines;

  // Setup internal variables
  int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                           
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                           
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                           
    Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;


  int N = num_splines;
#if defined HAVE_SSE2
  // We must pad to keep data align for SSE operations
  if (N & 1)
    N++;
#endif

  spline->x_stride = Ny*Nz*N;
  spline->y_stride = Nz*N;
  spline->z_stride = N;
  
#ifdef HAVE_POSIX_MEMALIGN
  posix_memalign ((void**)&spline->coefs, 64, ((size_t)sizeof(double)*Nx*Ny*Nz*N));
#else
  spline->coefs      = malloc ((size_t)sizeof(double)*Nx*Ny*Nz*N);
#endif

  spline->coefs_size=(size_t)Nx*(size_t)Ny*(size_t)Nz*(size_t)N;

#ifdef HAVE_SSE2
  init_sse_data();
#endif
  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_3d_d.\n");
    abort();
  }

  return spline;
}

void
set_multi_UBspline_3d_d (multi_UBspline_3d_d* spline, int num, double *data)
{
  int Mx = spline->x_grid.num;  
  int My = spline->y_grid.num; 
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                                   
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                                   
    Ny = My+2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                                   
    Nz = Mz+2;

  double *coefs = spline->coefs + num;
  intptr_t zs = spline->z_stride;

  // First, solve in the X-direction 
#pragma omp parallel for
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = (iy*Nz+iz)*zs;
      find_coefs_1d_d (spline->x_grid, spline->xBC, 
		       data+doffset,  (intptr_t)My*Mz,
		       coefs+coffset, (intptr_t)Ny*Nz*zs);
    }
  
  // Now, solve in the Y-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      intptr_t doffset = (ix*Ny*Nz + iz)*zs;
      intptr_t coffset = (ix*Ny*Nz + iz)*zs;
      find_coefs_1d_d (spline->y_grid, spline->yBC, 
		       coefs+doffset, (intptr_t)Nz*zs, 
		       coefs+coffset, (intptr_t)Nz*zs);
    }

  // Now, solve in the Z-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz*zs;
      intptr_t coffset = (ix*Ny+iy)*Nz*zs;
      find_coefs_1d_d (spline->z_grid, spline->zBC, 
		       coefs+doffset, (intptr_t)zs, 
		       coefs+coffset, (intptr_t)zs);
    }
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////    Double-Precision, Complex Creation Routines     ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions


void
destroy_multi_UBspline (Bspline *spline)
{
  free (spline->coefs);
  free (spline);
}
