//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include "blip_create.h"
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "config.h"
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif

#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <math.h>
#include <aligned_alloc.h>

void init_sse_data();

inline 
void* FFTAlign (void* ptr)
{
  size_t offset = 16 - (size_t)((size_t)ptr)&0x0f;
  return (void*) ((size_t)ptr+offset);
}

inline double dot (double a[3], double b[3])
{
  return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

// This function creates a single-precision real blip function from a
// set of plane-wave coefficients.  lattice is a 3x3 array specifying
// the lattice vectors.  The first lattice vector is given
// contiguously at latice[0], the second at lattice[3], and the third
// at lattice[6].  The next is a list of 3D G-vectors in the format:
// G_x[0] G_y[0] G_z[0], G_x[1], G_y[1], G_z[1],...
// Next, complex plane-wave coefficents are given, one for each
// G-vector.  Next, the number of G-vectors is given, followed by
// a factor which increases the density of the real-space grid.  A
// factor of 1.0 uses the minimum density to avoid aliasing.  Finally,
// the last parameter specifies whether to take the real or imaginary part.
// The spline is constructed to have domain [0,1) for x, y, and z coordinates. 
UBspline_3d_s*
create_blip_3d_s (double *lattice, double *Gvecs, 
		  complex_float *coefs, int numG,
		  double factor, bool useReal)
{
  int max_ix=0, max_iy=0, max_iz=0;
  int Mx, My, Mz;
  double twoPiInv = 1.0/(2.0*M_PI);
  for (int i=0; i<numG; i++) {
    double *G = Gvecs+3*i;
    int ix = round (twoPiInv * dot (lattice+0, G));
    int iy = round (twoPiInv * dot (lattice+3, G));
    int iz = round (twoPiInv * dot (lattice+6, G));
    if (abs(ix) > max_ix)   max_ix = ix;
    if (abs(iy) > max_iy)   max_iy = iy;
    if (abs(iz) > max_iz)   max_iz = iz;
  }
  Mx = 4*max_ix + 1;
  My = 4*max_iy + 1;
  Mz = 4*max_iz + 1;
  Mx = (int) ceil(factor*Mx);
  My = (int) ceil(factor*My);
  Mz = (int) ceil(factor*Mz);

  // FFTs are a little faster with even dimensions.
  if ((Mx%2)==1) Mx++;
  if ((My%2)==1) My++;
  if ((Mz%2)==1) Mz++;

  fprintf (stderr, "(Mx, My, Mz) = (%d, %d, %d)\n", Mx, My, Mz);

  // Now allocate space for FFT box
  complex_float *fft_box, *alloc_ptr;
  fft_box = aligned_alloc (sizeof(complex_float)*Mx*My*Mz, 16);

  // Create FFTW plan
  fftwf_plan plan = 
  fftwf_plan_dft_3d (Mx, My, Mz, (fftwf_complex*)fft_box, (fftwf_complex*)fft_box, 1,
		     FFTW_ESTIMATE);
  
  // Zero-out fft-box
  for (int i=0; i<Mx*My*Mz; i++)
    fft_box[i] = (complex_float)0.0f;
  
  // Now fill in fft box with coefficients in the right places
  double MxInv = 1.0/(double)Mx;
  double MyInv = 1.0/(double)My;
  double MzInv = 1.0/(double)Mz;
  double scale = 1.0/3.375;
  for (int i=0; i<numG; i++) {
    double *g = Gvecs+3*i;
    double G[3];
    G[0] = MxInv*(lattice[0]*g[0] + lattice[3]*g[1] + lattice[6]*g[2]);
    G[1] = MyInv*(lattice[1]*g[0] + lattice[4]*g[1] + lattice[7]*g[2]);
    G[2] = MzInv*(lattice[2]*g[0] + lattice[5]*g[1] + lattice[8]*g[2]);
    int ix = round (twoPiInv * dot (lattice+0, g));
    int iy = round (twoPiInv * dot (lattice+3, g));
    int iz = round (twoPiInv * dot (lattice+6, g));
    ix = (ix + Mx)%Mx;
    iy = (iy + My)%My;
    iz = (iz + Mz)%Mz;
    double gamma = 1.0;
    if (fabs(G[0]) > 1.0e-10)
      gamma *= (3.0/(G[0]*G[0]*G[0]*G[0])*(3.0 - 4.0*cos(G[0]) + cos(2.0*G[0])));
    else
      gamma *= 1.5;
    if (fabs(G[1]) > 1.0e-10)
      gamma *= (3.0/(G[1]*G[1]*G[1]*G[1])*(3.0 - 4.0*cos(G[1]) + cos(2.0*G[1])));
    else
      gamma *= 1.5;
    if (fabs(G[2]) > 1.0e-10)
      gamma *= (3.0/(G[2]*G[2]*G[2]*G[2])*(3.0 - 4.0*cos(G[2]) + cos(2.0*G[2])));
    else
      gamma *= 1.5;
    gamma *= scale;
    fft_box[(ix*My+iy)*Mz+iz] = coefs[i]/gamma;
  }
  
  // Execute the FFTW plan
  fftwf_execute (plan);
  // Destroy plan
  fftwf_destroy_plan (plan);

  // Now we have the coefficients in the FFT box.  We must allocate a
  // little bit larger box to hold the B-spline coefficients
  UBspline_3d_s* restrict spline = malloc (sizeof (UBspline_3d_s));
  spline->spcode = U3D;
  spline->tcode  = SINGLE_REAL;
  Ugrid x_grid, y_grid, z_grid;
  int Nx = Mx + 3;
  int Ny = My + 3;
  int Nz = Mz + 3;
  x_grid.start = 0.0;  x_grid.end = 1.0;  x_grid.num = Mx;  
  x_grid.delta = 1.0/(double)Mx;    x_grid.delta_inv = 1.0/x_grid.delta;
  y_grid.start = 0.0;  y_grid.end = 1.0;  y_grid.num = My;  
  y_grid.delta = 1.0/(double)My;    y_grid.delta_inv = 1.0/y_grid.delta;
  z_grid.start = 0.0;  z_grid.end = 1.0;  z_grid.num = Mz;  
  z_grid.delta = 1.0/(double)Mz;      z_grid.delta_inv = 1.0/z_grid.delta;
  spline->x_grid = x_grid;
  spline->y_grid = y_grid;
  spline->z_grid = z_grid;
  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;
  spline->xBC.lCode = PERIODIC;  spline->xBC.rCode = PERIODIC;
  spline->yBC.lCode = PERIODIC;  spline->yBC.rCode = PERIODIC;
  spline->zBC.lCode = PERIODIC;  spline->zBC.rCode = PERIODIC;
  
#ifndef HAVE_SSE2
  spline->coefs      = malloc (sizeof(float)*Nx*Ny*Nz);
#else
  posix_memalign ((void**)&spline->coefs, 16, sizeof(float)*Nx*Ny*Nz);
#endif 

  // Now copy data into spline coefficients, observing periodic boundary conditions
  for (int ix=0; ix<Nx; ix++) {
    int jx = (ix-1 + Mx)%Mx;
    for (int iy=0; iy < Ny; iy++) {
      int jy = (iy-1 + My)%My;
      for (int iz=0; iz < Nz; iz++) {
	int jz = (iz-1 + Mz)%Mz;
	if (useReal)
	  spline->coefs[(ix*Ny+iy)*Nz+iz] = 
	    crealf (fft_box[(jx*My+jy)*Mz+jz]);
	else
	  spline->coefs[(ix*Ny+iy)*Nz+iz] = 
	    cimagf (fft_box[(jx*My+jy)*Mz+jz]);	
      }
    }
  }
      
  //free (alloc_ptr);
  aligned_free (fft_box);

  init_sse_data();
  return spline;
}
