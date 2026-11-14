/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include "multi_bspline_create.h"
#include <cstdlib>
#include <cstdio>
#include <cstdint>

void init_sse_data();

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
  multi_UBspline_1d_s* restrict spline = static_cast<multi_UBspline_1d_s*>(std::malloc (sizeof(multi_UBspline_1d_s)));
  if (!spline) {
    std::fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_1d_s.\n");
    std::abort();
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
  spline->coefs = static_cast<float*>(std::malloc (sizeof(float)*Nx*N));
#else
  void* coefs_ptr = nullptr;
  if (posix_memalign (&coefs_ptr, 64, (sizeof(float)*Nx*N)) != 0)
    coefs_ptr = nullptr;
  spline->coefs = static_cast<float*>(coefs_ptr);
#endif
  spline->coefs_size=static_cast<size_t>(Nx)*static_cast<size_t>(N);
#ifdef HAVE_SSE
  init_sse_data();    
#endif
  if (!spline->coefs) {
    std::fprintf (stderr, "Out of memory allocating spline coefficient in create_multi_UBspline_1d_s.\n");
    std::abort();
  }


  return spline;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////     Double-Precision, Real Creation Routines       ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

multi_UBspline_1d_d*
create_multi_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, int num_splines)
{
  // Create new spline
  multi_UBspline_1d_d* restrict spline = static_cast<multi_UBspline_1d_d*>(std::malloc (sizeof(multi_UBspline_1d_d)));
  if (!spline) {
    std::fprintf (stderr, "Out of memory allocating spline in create_multi_UBspline_1d_d.\n");
    std::abort();
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
  spline->coefs = static_cast<double*>(std::malloc (sizeof(double)*Nx*N));
#else
  void* coefs_ptr = nullptr;
  if (posix_memalign (&coefs_ptr, 64, sizeof(double)*Nx*N) != 0)
    coefs_ptr = nullptr;
  spline->coefs = static_cast<double*>(coefs_ptr);
#endif
  spline->coefs_size=static_cast<size_t>(Nx)*static_cast<size_t>(N);
#ifdef HAVE_SSE2
  init_sse_data();
#endif
  if (!spline->coefs) {
    std::fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_1d_d.\n");
    std::abort();
  }

  return spline;
}

void
destroy_multi_UBspline (Bspline *spline)
{
  std::free (spline->coefs);
  std::free (spline);
}
