//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#include<cassert>
#include <complex>
#include<cuda.h>
#include <thrust/complex.h>
#include<cuda_runtime.h>
#include "AFQMC/Kernels/cuda_settings.h"
#define QMC_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{

// very sloppy, needs improvement!!!!
template<typename T, typename T2>
__global__ void kernel_vKKwij_to_vwKiKj( int nwalk, int nkpts, int nmo_max, int nmo_tot, 
                                        bool* kk, int* nmo, int* nmo0, T const* A, T2 * B)
{
  int Ki = blockIdx.x;  
  int Kj = blockIdx.y;  
  int nw = blockIdx.z;
  if(Ki >= nkpts || Kj >= nkpts || nw >= nwalk) return;  
  int ni0 = nmo0[Ki]; 
  int nj0 = nmo0[Kj]; 
  int ni = nmo[Ki]; 
  int nj = nmo[Kj];
  
  T2* B_(B + nw*nmo_tot*nmo_tot + ni0*nmo_tot + nj0); 
  T const* A_(A + ((Ki*nkpts + Kj)*nwalk + nw)*nmo_max*nmo_max);

  if(threadIdx.x >= ni) return; 
  if(threadIdx.y >= nj) return; 
  
  if(kk[Ki*nkpts+Kj]) { // transpose
    for(int i=threadIdx.x; i<ni; i+=blockDim.x) { 
      for(int j=threadIdx.y; j<nj; j+=blockDim.y) {
        B_[ i*nmo_tot+j ] = static_cast<T2>(A_[ j*nmo_max+i ]);
      }
    }
  } else { // copy
    for(int i=threadIdx.x; i<ni; i+=blockDim.x) { 
      for(int j=threadIdx.y; j<nj; j+=blockDim.y) {
        B_[ i*nmo_tot+j ] = static_cast<T2>(A_[ i*nmo_max+j ]);
      }
    }
  }  
}

template<typename T, typename T2>
__global__ void kernel_vKKwij_to_vwKiKj( int nwalk, int nkpts, int nmo_max, int nmo_tot, 
                                        bool* kk, int* nmo, int* nmo0, thrust::complex<T> const* A, thrust::complex<T2> * B)
{
// use shared memory for transpose
//  __shared__ thrust::complex<T2> tile[TILE_DIM][TILE_DIM+1];
  int Ki = blockIdx.x;
  int Kj = blockIdx.y;
  int nw = blockIdx.z;
  if(Ki >= nkpts || Kj >= nkpts || nw >= nwalk) return;
  int ni0 = nmo0[Ki];
  int nj0 = nmo0[Kj];
  int ni = nmo[Ki];  
  int nj = nmo[Kj];  
            
  thrust::complex<T2>* B_(B + nw*nmo_tot*nmo_tot + ni0*nmo_tot + nj0);  
  thrust::complex<T> const* A_(A + ((Ki*nkpts + Kj)*nwalk + nw)*nmo_max*nmo_max);

  if(threadIdx.x >= ni) return; 
  if(threadIdx.y >= nj) return; 

  if(kk[Ki*nkpts+Kj]) { // transpose
    for(int i=threadIdx.x; i<ni; i+=blockDim.x) { 
      for(int j=threadIdx.y; j<nj; j+=blockDim.y) {
        B_[ i*nmo_tot+j ] = static_cast<thrust::complex<T2>>(A_[ j*nmo_max+i ]);
      }
    }
  } else {  // copy
    for(int i=threadIdx.x; i<ni; i+=blockDim.x) { 
      for(int j=threadIdx.y; j<nj; j+=blockDim.y) {
        B_[ i*nmo_tot+j ] = static_cast<thrust::complex<T2>>(A_[ i*nmo_max+j ]);
      }
    }
  }
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                               bool* kk,int* nmo, int* nmo0, double const* A, double * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,A,B);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                               bool* kk, int* nmo, int* nmo0, float const* A, float * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,A,B);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                               bool* kk, int* nmo, int* nmo0, float const* A, double * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,A,B);
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                                bool* kk, int* nmo, int* nmo0, std::complex<double> const* A, std::complex<double> * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,
                reinterpret_cast<thrust::complex<double> const*>(A),
                reinterpret_cast<thrust::complex<double> *>(B));
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                                bool* kk, int* nmo, int* nmo0, std::complex<float> const* A, std::complex<float> * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,
                reinterpret_cast<thrust::complex<float> const*>(A),
                reinterpret_cast<thrust::complex<float> *>(B));
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void vKKwij_to_vwKiKj(int nwalk, int nkpts, int nmo_max, int nmo_tot,
                                bool* kk, int* nmo, int* nmo0, std::complex<float> const* A, std::complex<double> * B)
{
  int xblock_dim = 8;
  int yblock_dim = 8;
  dim3 block_dim(xblock_dim,yblock_dim,1);
  dim3 grid_dim(nkpts,nkpts,nwalk);
  kernel_vKKwij_to_vwKiKj<<<grid_dim, block_dim>>>(nwalk,nkpts,nmo_max,nmo_tot,kk,nmo,nmo0,
                reinterpret_cast<thrust::complex<float> const*>(A),
                reinterpret_cast<thrust::complex<double> *>(B));
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


}
