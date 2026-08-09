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
destroy_multi_UBspline (Bspline *spline)
{
  free (spline->coefs);
  free (spline);
}
