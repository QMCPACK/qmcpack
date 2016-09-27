//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>

#include "bspline_base.h"
#include "bspline_structs.h"
#include "bspline_structs_cuda.h"

__device__ double Bcuda[48];
__constant__ float  Acuda[48];

// #include "bspline_cuda_s_impl.h"
// #include "bspline_cuda_c_impl.h"
// #include "bspline_cuda_d_impl.h"
// #include "bspline_cuda_z_impl.h"

extern "C" UBspline_3d_c_cuda*
create_UBspline_3d_c_cuda (UBspline_3d_c* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  UBspline_3d_c_cuda *cuda_spline =
    (UBspline_3d_c_cuda*) malloc (sizeof (UBspline_3d_c_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32)*32;

  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(std::complex<float>);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  std::complex<float> *spline_buff = (std::complex<float>*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++) {
      for (int iz=0; iz<Nz; iz++) 
	spline_buff[ix*cuda_spline->stride.x +
		    iy*cuda_spline->stride.y + iz] =
	  spline->coefs[ix*spline->x_stride +
			iy*spline->y_stride +iz]; 
      for (int isp=Nz; isp < N; isp++) {
	spline_buff[ix*cuda_spline->stride.x +
		    iy*cuda_spline->stride.y + isp] = 0.0;
      }
    }
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  free(spline_buff);
  
  cuda_spline->stride.x = 2*Ny*N;
  cuda_spline->stride.y = 2*N;

  return cuda_spline;
}


extern "C" UBspline_3d_c_cuda*
create_UBspline_3d_c_cuda_conv (UBspline_3d_z* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  UBspline_3d_c_cuda *cuda_spline =
    (UBspline_3d_c_cuda*) malloc (sizeof (UBspline_3d_c_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32) * 32;
  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(std::complex<float>);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  std::complex<float> *spline_buff = (std::complex<float>*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) {
	std::complex<double> z = spline->coefs[ix*spline->x_stride +
					       iy*spline->y_stride + iz];
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y + iz] = std::complex<float>(z.real(), z.imag());
	for (int iz=Nz; iz < N; iz++) 
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y + iz] = 0.0;
      }

	

  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  cuda_spline->stride.x = 2*Ny*N;
  cuda_spline->stride.y = 2*N;

  return cuda_spline;
}




extern "C" UBspline_3d_s_cuda*
create_UBspline_3d_s_cuda (UBspline_3d_s* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  UBspline_3d_s_cuda *cuda_spline =
    (UBspline_3d_s_cuda*) malloc (sizeof (UBspline_3d_s_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32)*32;

  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;
  cuda_spline->stride.z = 1;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	spline_buff[ix*cuda_spline->stride.x +
		    iy*cuda_spline->stride.y + iz] = 
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride + iz];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}



extern "C" UBspline_3d_s_cuda*
create_UBspline_3d_s_cuda_conv (UBspline_3d_d* spline)
{
  float A_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Acuda, A_h, 48*sizeof(float), 0, cudaMemcpyHostToDevice);

  UBspline_3d_s_cuda *cuda_spline =
    (UBspline_3d_s_cuda*) malloc (sizeof (UBspline_3d_s_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32)*32;
  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;
  cuda_spline->stride.z = 1;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to allocate %ld memory for GPU spline coefficients.  Error %s\n",
	     size, cudaGetErrorString(err));
    abort();
  }
  
  float *spline_buff = (float*)malloc(size);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y + iz] = 
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride + iz];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to copy spline to GPU memory.  Error:  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
  free(spline_buff);

  return cuda_spline;
}




extern "C" UBspline_3d_d_cuda*
create_UBspline_3d_d_cuda (UBspline_3d_d* spline)
{
  double B_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Bcuda, B_h, 48*sizeof(double), 0, cudaMemcpyHostToDevice);

  UBspline_3d_d_cuda *cuda_spline =
    (UBspline_3d_d_cuda*) malloc (sizeof (UBspline_3d_d_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32)*32;
  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;
  cuda_spline->stride.z = 1;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(double);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  double *spline_buff = (double*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y + iz] = 
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride + iz];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}



extern "C" UBspline_3d_z_cuda*
create_UBspline_3d_z_cuda (UBspline_3d_z* spline)
{
  double B_h[48] = { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
		     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
		    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
		     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0,
		         0.0,     -0.5,      1.0,    -0.5,
		         0.0,      1.5,     -2.0,     0.0,
		         0.0,     -1.5,      1.0,     0.5,
		         0.0,      0.5,      0.0,     0.0,
		         0.0,      0.0,     -1.0,     1.0,
		         0.0,      0.0,      3.0,    -2.0,
		         0.0,      0.0,     -3.0,     1.0,
		         0.0,      0.0,      1.0,     0.0 };

  cudaMemcpyToSymbol(Bcuda, B_h, 48*sizeof(double), 0, cudaMemcpyHostToDevice);

  UBspline_3d_z_cuda *cuda_spline =
    (UBspline_3d_z_cuda*) malloc (sizeof (UBspline_3d_z_cuda));
  
  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = ((Nz+31)/32)*32;
  cuda_spline->stride.x = Ny*N;
  cuda_spline->stride.y = N;  
  cuda_spline->stride.z = 1;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*N*sizeof(std::complex<double>);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  std::complex<double> *spline_buff = (std::complex<double>*)malloc(size);

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y + iz] = 
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride + iz];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  cuda_spline->stride.x = 2*Ny*N;
  cuda_spline->stride.y = 2*N;
  cuda_spline->stride.z = 2;

  free(spline_buff);

  return cuda_spline;
}
