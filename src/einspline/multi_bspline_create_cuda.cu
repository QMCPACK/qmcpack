//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//		      Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>

#include "multi_bspline.h"
#include "multi_bspline_structs_cuda.h"

__device__ double Bcuda[48];
__constant__ float  Acuda[48];

#include "multi_bspline_cuda_s_impl.h"
#include "multi_bspline_cuda_c_impl.h"
#include "multi_bspline_cuda_d_impl.h"
#include "multi_bspline_cuda_z_impl.h"

#include <CUDA/gpu_vector.h>

#define COALLESCED_SIZE 16

#define UnifiedVirtualAddressing

extern "C" multi_UBspline_1d_s_cuda*
create_multi_UBspline_1d_s_cuda (multi_UBspline_1d_s* spline)
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

  multi_UBspline_1d_s_cuda *cuda_spline =
    (multi_UBspline_1d_s_cuda*) malloc (sizeof (multi_UBspline_1d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int N = spline->num_splines;

  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride = N;
  cuda_spline->gridInv = spline->x_grid.delta_inv;
  cuda_spline->dim = spline->x_grid.num;

  size_t size = Nx*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }

  for (int ix=0; ix<Nx; ix++)
    for (int isp=0; isp<spline->num_splines; isp++) 
      spline_buff[ix*cuda_spline->stride + isp] =
	spline->coefs[ix*spline->x_stride + isp];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  
  free(spline_buff);
  
  return cuda_spline;
}


extern "C" multi_UBspline_1d_s_cuda*
create_multi_UBspline_1d_s_cuda_conv (multi_UBspline_1d_d* spline)
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

  multi_UBspline_1d_s_cuda *cuda_spline =
    (multi_UBspline_1d_s_cuda*) malloc (sizeof (multi_UBspline_1d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int N = spline->num_splines;

  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride = N;
  cuda_spline->gridInv = spline->x_grid.delta_inv;
  cuda_spline->dim = spline->x_grid.num;

  size_t size = Nx*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }


  for (int ix=0; ix<Nx; ix++)
    for (int isp=0; isp<spline->num_splines; isp++) 
      spline_buff[ix*cuda_spline->stride + isp] = 
	(float)spline->coefs[ix*spline->x_stride + isp];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  
  free(spline_buff);
  
  return cuda_spline;
}



extern "C" multi_UBspline_1d_c_cuda*
create_multi_UBspline_1d_c_cuda (multi_UBspline_1d_c* spline)
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

  multi_UBspline_1d_c_cuda *cuda_spline =
    (multi_UBspline_1d_c_cuda*) malloc (sizeof (multi_UBspline_1d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int N = spline->num_splines;

  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride = N;
  cuda_spline->gridInv = spline->x_grid.delta_inv;
  cuda_spline->dim = spline->x_grid.num;

  size_t size = Nx*N*sizeof(complex_float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  complex_float *spline_buff = (complex_float*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }


  for (int ix=0; ix<Nx; ix++)
    for (int isp=0; isp<spline->num_splines; isp++) 
      spline_buff[ix*cuda_spline->stride + isp] =
	spline->coefs[ix*spline->x_stride + isp];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  
  free(spline_buff);
  
  return cuda_spline;
}


extern "C" multi_UBspline_1d_c_cuda*
create_multi_UBspline_1d_c_cuda_conv (multi_UBspline_1d_z* spline)
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
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Error copying A matrix to GPU constant memory:  Erorr = %s\n",
	     cudaGetErrorString(err));
    abort();
  }

  multi_UBspline_1d_c_cuda *cuda_spline =
    (multi_UBspline_1d_c_cuda*) malloc (sizeof (multi_UBspline_1d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int N = spline->num_splines;

  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride = N;
  cuda_spline->gridInv = spline->x_grid.delta_inv;
  cuda_spline->dim = spline->x_grid.num;

  size_t size = Nx*N*sizeof(complex_float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  complex_float *spline_buff = (complex_float*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }


  for (int ix=0; ix<Nx; ix++)
    for (int isp=0; isp<spline->num_splines; isp++) 
      spline_buff[ix*cuda_spline->stride + isp] =
	spline->coefs[ix*spline->x_stride + isp];
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  
  free(spline_buff);
  
  return cuda_spline;
}




extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda (multi_UBspline_3d_c* spline)
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

  multi_UBspline_3d_c_cuda *cuda_spline =
    (multi_UBspline_3d_c_cuda*) malloc (sizeof (multi_UBspline_3d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  cuda_spline->host_Nx_offset = spline->x_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(std::complex<float>);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  std::complex<float> *spline_buff = (std::complex<float>*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }


  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) {
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
	for (int isp=spline->num_splines; isp < N; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] = 0.0;
	}

      }
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  free(spline_buff);

  cuda_spline->stride.x = 2*Ny*Nz*N;
  cuda_spline->stride.y = 2*Nz*N;
  cuda_spline->stride.z = 2*N;


  return cuda_spline;
}


extern "C" multi_UBspline_3d_c_cuda*
create_multi_UBspline_3d_c_cuda_conv (multi_UBspline_3d_z* spline)
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

  multi_UBspline_3d_c_cuda *cuda_spline =
    (multi_UBspline_3d_c_cuda*) malloc (sizeof (multi_UBspline_3d_c_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  cuda_spline->host_Nx_offset = spline->x_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(std::complex<float>);

#ifdef UnifiedVirtualAddressing

  size_t max_GPU_spine_size = gpu::MaxGPUSpineSizeMB;
  if (!gpu::rank && size > (max_GPU_spine_size<<20)) {
    fprintf (stderr, "Required spline table %ld MB on the GPU is larger than the limit %ld MB.\n", size>>20, max_GPU_spine_size);
    fprintf (stderr, "Partial table is stored on the host memory.\n");
  }
  int Nx_host;
  int Nx_GPU = (max_GPU_spine_size<<20)/(cuda_spline->stride.x*sizeof(std::complex<float>));
  if (Nx_GPU >= Nx)
    { Nx_GPU = Nx; Nx_host = 0; cuda_spline->host_Nx_offset = Nx_GPU - 3; }
  else if (Nx_GPU < 4)
    { Nx_GPU = 0; Nx_host = Nx; cuda_spline->host_Nx_offset = 0; } 
  else
    { Nx_host = Nx - Nx_GPU + 3; cuda_spline->host_Nx_offset = Nx_GPU - 3; }
  if(!gpu::rank) fprintf (stderr, "Using Nx %d on GPU %d, on host %d, host_Nx_offset %d\n", Nx, Nx_GPU, Nx_host, cuda_spline->host_Nx_offset);

  size_t size_GPU = Nx_GPU*cuda_spline->stride.x*sizeof(std::complex<float>);
  size_t size_host = Nx_host*cuda_spline->stride.x*sizeof(std::complex<float>);

  gpu::MaxGPUSpineSizeMB = gpu::MaxGPUSpineSizeMB - ((size_t)Nx_GPU)*cuda_spline->stride.x*sizeof(std::complex<float>)/1024/1024;

  if (Nx_GPU)
  {
    cuda_spline->coefs = (complex_float *) gpu::cuda_memory_manager.allocate(size_GPU, "SPO_multi_UBspline_3d_c_cuda");
    std::complex<float> *spline_buff = (std::complex<float>*)malloc(size_GPU);
    if (!spline_buff) {
      fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
      abort();
    }

    for (int ix=0; ix<Nx_GPU; ix++)
      for (int iy=0; iy<Ny; iy++)
        for (int iz=0; iz<Nz; iz++) {
          for (int isp=0; isp<spline->num_splines; isp++) {
            std::complex<double> z = spline->coefs[ix*spline->x_stride + iy*spline->y_stride + iz*spline->z_stride + isp];
            spline_buff[ix*cuda_spline->stride.x + iy*cuda_spline->stride.y + iz*cuda_spline->stride.z + isp] 
              = std::complex<float>(z.real(), z.imag());
          }
          for (int isp=spline->num_splines; isp < N; isp++) 
            spline_buff[ix*cuda_spline->stride.x + iy*cuda_spline->stride.y + iz*cuda_spline->stride.z + isp] = 0.0;
        }

    cudaMemcpy(cuda_spline->coefs, spline_buff, size_GPU, cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf (stderr, "Failed to copy spline to GPU memory.  Error:  %s\n", cudaGetErrorString(err));
      abort();
    }
    free(spline_buff);
  }

  if (Nx_host)
  {
    std::complex<float> *spline_buff;
    cudaHostAlloc(&spline_buff, size_host, cudaHostAllocMapped);
    if (!spline_buff) {
      fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
      abort();
    }

    for (int ix=0; ix<Nx_host; ix++)
      for (int iy=0; iy<Ny; iy++)
        for (int iz=0; iz<Nz; iz++) {
          for (int isp=0; isp<spline->num_splines; isp++) {
  	      std::complex<double> z = spline->coefs[(ix+Nx-Nx_host)*spline->x_stride + iy*spline->y_stride + iz*spline->z_stride + isp];
  	      spline_buff[ix*cuda_spline->stride.x + iy*cuda_spline->stride.y + iz*cuda_spline->stride.z + isp] 
              = std::complex<float>(z.real(), z.imag());
  	    }
          for (int isp=spline->num_splines; isp < N; isp++) 
  	      spline_buff[ix*cuda_spline->stride.x + iy*cuda_spline->stride.y + iz*cuda_spline->stride.z + isp] = 0.0;
          }
    cuda_spline->coefs_host = spline_buff;
  }
#else
  cuda_spline->coefs = (complex_float *) gpu::cuda_memory_manager.allocate(size, "SPO_multi_UBspline_3d_c_cuda");
  /*
  cudaMalloc((void**)&(cuda_spline->coefs), size);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to allocate %ld memory for GPU spline coefficients.  Error %s\n",
	     size, cudaGetErrorString(err));
    abort();
  }*/
  //fprintf (stdout, "Succeeded to allocate %ld bytes memory for GPU spline coefficients.\n", size);
 
  std::complex<float> *spline_buff = (std::complex<float>*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) {
	for (int isp=0; isp<spline->num_splines; isp++) {
	  std::complex<double> z = spline->coefs[ix*spline->x_stride +
						 iy*spline->y_stride +
						 iz*spline->z_stride + isp];
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] = std::complex<float>(z.real(), z.imag());
	}
	for (int isp=spline->num_splines; isp < N; isp++) 
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] = 0.0;
      }

  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to copy spline to GPU memory.  Error:  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
  free(spline_buff);
#endif

  cuda_spline->stride.x = 2*Ny*Nz*N;
  cuda_spline->stride.y = 2*Nz*N;
  cuda_spline->stride.z = 2*N;


  return cuda_spline;
}




extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda (multi_UBspline_3d_s* spline)
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

  multi_UBspline_3d_s_cuda *cuda_spline =
    (multi_UBspline_3d_s_cuda*) malloc (sizeof (multi_UBspline_3d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cudaMalloc((void**)&(cuda_spline->coefs), size);
  
  float *spline_buff = (float*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }


  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}



extern "C" multi_UBspline_3d_s_cuda*
create_multi_UBspline_3d_s_cuda_conv (multi_UBspline_3d_d* spline)
{
  // fprintf (stderr, "In create_multi_UBspline_3d_s_cuda_conv.\n");
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

  multi_UBspline_3d_s_cuda *cuda_spline =
    (multi_UBspline_3d_s_cuda*) malloc (sizeof (multi_UBspline_3d_s_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(float);

  cuda_spline->coefs = (float *) gpu::cuda_memory_manager.allocate(size, "SPO_multi_UBspline_3d_s_cuda");
  float *spline_buff = (float *)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	  // if (isnan (spline->coefs[ix*spline->x_stride +
	  // 			   iy*spline->y_stride +
	  // 			   iz*spline->z_stride + isp]))
	  //    fprintf (stderr, "NAN at ix=%d iy=%d iz=%d isp=%d\n",
	  //    	     ix,iy,iz,isp);
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to copy spline to GPU memory.  Error:  %s\n",
	     cudaGetErrorString(err));
    abort();
  }
  free(spline_buff);

  return cuda_spline;
}




extern "C" multi_UBspline_3d_d_cuda*
create_multi_UBspline_3d_d_cuda (multi_UBspline_3d_d* spline)
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

  multi_UBspline_3d_d_cuda *cuda_spline =
    (multi_UBspline_3d_d_cuda*) malloc (sizeof (multi_UBspline_3d_d_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(double);

  cuda_spline->coefs = (double *) gpu::cuda_memory_manager.allocate(size, "SPO_multi_UBspline_3d_d_cuda");
  double *spline_buff = (double *)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  free(spline_buff);

  return cuda_spline;
}



extern "C" multi_UBspline_3d_z_cuda*
create_multi_UBspline_3d_z_cuda (multi_UBspline_3d_z* spline)
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

  multi_UBspline_3d_z_cuda *cuda_spline =
    (multi_UBspline_3d_z_cuda*) malloc (sizeof (multi_UBspline_3d_z_cuda));
  
  cuda_spline->num_splines = spline->num_splines;

  int Nx = spline->x_grid.num+3;
  int Ny = spline->y_grid.num+3;
  int Nz = spline->z_grid.num+3;

  int N = spline->num_splines;
  if ((N%COALLESCED_SIZE) != 0)
    N += COALLESCED_SIZE - (N%COALLESCED_SIZE);
  cuda_spline->stride.x = Ny*Nz*N;
  cuda_spline->stride.y = Nz*N;
  cuda_spline->stride.z = N;

  cuda_spline->gridInv.x = spline->x_grid.delta_inv;
  cuda_spline->gridInv.y = spline->y_grid.delta_inv;
  cuda_spline->gridInv.z = spline->z_grid.delta_inv;

  cuda_spline->dim.x = spline->x_grid.num;
  cuda_spline->dim.y = spline->y_grid.num;
  cuda_spline->dim.z = spline->z_grid.num;

  size_t size = Nx*Ny*Nz*N*sizeof(std::complex<double>);

  cuda_spline->coefs = (complex_double *) gpu::cuda_memory_manager.allocate(size, "SPO_multi_UBspline_3d_z_cuda");
  /*
  cudaMalloc((void**)&(cuda_spline->coefs), size);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "Failed to allocate %ld memory for GPU spline coefficients.  Error %s\n",
	     size, cudaGetErrorString(err));
    abort();
  }*/
  
  std::complex<double> *spline_buff = (std::complex<double>*)malloc(size);
  if (!spline_buff) {
    fprintf (stderr, "Failed to allocate memory for temporary spline buffer.\n");
    abort();
  }

  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++) 
	for (int isp=0; isp<spline->num_splines; isp++) {
	  spline_buff[ix*cuda_spline->stride.x +
		      iy*cuda_spline->stride.y +
		      iz*cuda_spline->stride.z + isp] =
	    spline->coefs[ix*spline->x_stride +
			  iy*spline->y_stride +
			  iz*spline->z_stride + isp];
	}
  cudaMemcpy(cuda_spline->coefs, spline_buff, size, cudaMemcpyHostToDevice);

  cuda_spline->stride.x = 2*Ny*Nz*N;
  cuda_spline->stride.y = 2*Nz*N;
  cuda_spline->stride.z = 2*N;

  free(spline_buff);

  return cuda_spline;
}
