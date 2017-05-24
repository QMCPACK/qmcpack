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
/** @file einspline_allocator.c
 *
 * Gather only 3d_d/3d_s allocation routines.
 * Original implementations are einspline/ *_create.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "config.h"
#include "einspline/bspline.h"
#include "spline2/einspline_allocator.h"

//forward declaration of the functions
void find_coefs_1d_d (Ugrid grid, BCtype_d bc, double *data,  intptr_t dstride,
    double *coefs, intptr_t cstride);

void solve_deriv_interp_1d_s (float bands[], float coefs[], int M, int
    cstride);

void solve_antiperiodic_interp_1d_s (float bands[], float coefs[], int M, int
    cstride);

void find_coefs_1d_s (Ugrid grid, BCtype_s bc, float *data,  intptr_t dstride,
    float *coefs, intptr_t cstride);


#if defined(__INTEL_COMPILER)

inline void *
einspline_alloc (size_t size, size_t alignment)
{
  return _mm_malloc(size,alignment);
}

inline void
einspline_free (void *ptr)
{
  _mm_free(ptr);
}

#elif defined(HAVE_POSIX_MEMALIGN)

int posix_memalign(void **memptr, size_t alignment, size_t size);

inline void *
einspline_alloc (size_t size, size_t alignment)
{
  void *ptr;
  posix_memalign (&ptr, alignment, size);
  return ptr;
}

inline void
einspline_free (void *ptr)
{
  free (ptr);
}

#else

inline void *
einspline_alloc (size_t size, size_t alignment)
{
  size += (alignment-1)+sizeof(void*);
  void *ptr = malloc (size);
  if (ptr == NULL)
    return NULL;
  else
  {
    void *shifted = (char *)ptr + sizeof(void*); // make room to save original pointer
    size_t offset = alignment - (size_t)shifted%(size_t)alignment;
    void *aligned = (char *)shifted + offset;
    *((void**)aligned-1) = ptr;
    return aligned;
  }
}

inline void
einspline_free (void *aligned)
{
  void *ptr = *((void**)aligned-1);
  free (ptr);
}
#endif


multi_UBspline_3d_s*
einspline_create_multi_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
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

  const int ND=QMC_CLINE/sizeof(float);
  const size_t N=((num_splines+ND-1)/ND)*ND;
  //int N= (num_splines%ND) ? (num_splines+ND-num_splines%ND) : num_splines;
  spline->x_stride = (size_t)Ny*(size_t)Nz*(size_t)N;
  spline->y_stride = (size_t)Nz*N;
  spline->z_stride = N;

  spline->coefs_size=(size_t)Nx*spline->x_stride;
  spline->coefs = (float*)einspline_alloc(sizeof(float)*spline->coefs_size,QMC_CLINE);
  //printf("Einepline allocator %d %d %d %zd (%d)  %zu %d\n",Nx,Ny,Nz,N,num_splines,spline->coefs_size,QMC_CLINE);

  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_3d_s.\n");
    abort();
  }

#if 0
  //test first-touch later
  const size_t xs = spline->x_stride;
  const size_t ys = spline->y_stride;
  const size_t zs = spline->z_stride;

  const float czero=0;
#pragma omp parallel for collapse(3)
  for(size_t i=0; i<Nx; ++i)
    for(size_t j=0; j<Ny; ++j)
      for(size_t k=0; k<Nz; ++k)
      {
        float* restrict coefs = spline->coefs + i*xs + j*ys + k*zs; 
        for(size_t s=0; s<N; ++s)
          coefs[s]=czero;
      }
#endif

  return spline;
}

multi_UBspline_3d_d*
einspline_create_multi_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
			    BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
			    int num_splines)
{
  // Create new spline
  multi_UBspline_3d_d* restrict spline = malloc (sizeof(multi_UBspline_3d_d));

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

  const int ND=QMC_CLINE/sizeof(double);
  int N= (num_splines%ND) ? (num_splines+ND-num_splines%ND) : num_splines;

  spline->x_stride = (size_t)Ny*(size_t)Nz*(size_t)N;
  spline->y_stride = Nz*N;
  spline->z_stride = N;

  spline->coefs_size=(size_t)Nx*spline->x_stride;
  spline->coefs=(double*)einspline_alloc(sizeof(double)*spline->coefs_size,QMC_CLINE);
  
  if (!spline->coefs) {
    fprintf (stderr, "Out of memory allocating spline coefficients in create_multi_UBspline_3d_d.\n");
    abort();
  }

  return spline;
}

UBspline_3d_d*
einspline_create_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
		      double *data)
{
  // Create new spline
  UBspline_3d_d* restrict spline = malloc (sizeof(UBspline_3d_d));
  spline->spcode = U3D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 

  // Setup internal variables
  int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs_size=(size_t)Nx*(size_t)Ny*(size_t)Nz;

  spline->coefs= (double*)einspline_alloc (sizeof(double)*spline->coefs_size,QMC_CLINE);

  if(data != NULL) // only data is provided
  {
  // First, solve in the X-direction 
#pragma omp parallel for
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_d (spline->x_grid, xBC, data+doffset, My*Mz,
		       spline->coefs+coffset, Ny*Nz);
    }
  
  // Now, solve in the Y-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_d (spline->y_grid, yBC, spline->coefs+doffset, Nz, 
		       spline->coefs+coffset, Nz);
    }

  // Now, solve in the Z-direction
#pragma omp parallel for
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_d (spline->z_grid, zBC, spline->coefs+doffset, 1, 
		       spline->coefs+coffset, 1);
    }
  }
  return spline;
}


UBspline_3d_s*
einspline_create_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
		      BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
		      float *data)
{
  // Create new spline
  UBspline_3d_s* spline = malloc (sizeof(UBspline_3d_s));
  spline->spcode = U3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  // Setup internal variables
  int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs_size=(size_t)Nx*(size_t)Ny*(size_t)Nz;
  spline->coefs= (float*)einspline_alloc (sizeof(float)*spline->coefs_size,QMC_CLINE);

  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) 
    for (int iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_s (spline->x_grid, xBC, data+doffset, My*Mz,
		       spline->coefs+coffset, Ny*Nz);
    }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_s (spline->y_grid, yBC, spline->coefs+doffset, Nz, 
		       spline->coefs+coffset, Nz);
    }

  // Now, solve in the Z-direction
  for (int ix=0; ix<Nx; ix++) 
    for (int iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_s (spline->z_grid, zBC, spline->coefs+doffset, 1, 
		       spline->coefs+coffset, 1);
    }
  return spline;
}
