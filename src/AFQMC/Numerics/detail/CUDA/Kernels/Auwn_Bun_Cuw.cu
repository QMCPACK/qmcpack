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
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels 
{

// C[u][w] = alpha * sum_a A[u][w][a] * B[u][a]
template<typename T>
__global__ void kernel_Auwn_Bun_Cuw( int nu, int nw, int na, 
            thrust::complex<T> const alpha, 
            thrust::complex<T> const* A, thrust::complex<T> const* B, 
            thrust::complex<T>* C)
{
    int nu_per_block = blockDim.x;
    int u = blockIdx.x*nu_per_block + threadIdx.x;
    int w = threadIdx.y; 

    if( (u<nu) && (w<nw) ) {
      thrust::complex<T> Cuw=0;
      thrust::complex<T> const* A_( A + (u*nw + w)*na );
      thrust::complex<T> const* B_( B + u*na );
      for(int a=0; a<na; ++a, ++A_, ++B_) Cuw += (*A_) * (*B_);
      C[ u*nw + w ] = alpha*Cuw;
    }
}

// NOT OPTIMAL: poor memory access pattern
template<typename T>
__global__ void kernel_Awiu_Biu_Cuw( int nu, int nw, int ni,
            thrust::complex<T> const alpha,
            thrust::complex<T> const* A, T const* B, int ldb,
            thrust::complex<T>* C, int ldc)
{
    int nu_per_block = blockDim.x;
    int u = blockIdx.x*nu_per_block + threadIdx.x;
    int w = threadIdx.y;

    if( (u<nu) && (w<nw) ) {
      thrust::complex<T> Cuw=0;
      thrust::complex<T> const* A_( A + w*ni*nu + u );
      T const* B_( B + u );
      for(int i=0; i<ni; ++i, A_+=nu, B_+=ldb) Cuw += (*A_) * (*B_);
      C[ u*ldc + w ] = alpha*Cuw;
    }
}

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
template<typename T>
__global__ void kernel_Awiu_Biu_Cuw( int nu, int nw, int ni,
            thrust::complex<T> const alpha,
            thrust::complex<T> const* A, thrust::complex<T> const* B, int ldb,
            thrust::complex<T>* C, int ldc)
{
    int nu_per_block = blockDim.x;
    int u = blockIdx.x*nu_per_block + threadIdx.x;
    int w = threadIdx.y;

    if( (u<nu) && (w<nw) ) {
      thrust::complex<T> Cuw=0;
      thrust::complex<T> const* A_( A + w*ni*nu + u );
      thrust::complex<T> const* B_( B + u );
      for(int i=0; i<ni; ++i, A_+=nu, B_+=ldb) Cuw += (*A_) * (*B_);
      C[ u*ldc + w ] = alpha*Cuw;
    }
}

// element-wise C[k][i][j] = A[i][j] * B[j][k]
template<typename T>
__global__ void kernel_element_wise_Aij_Bjk_Ckij( char transA, int ni, int nj, int nk,
            T const* A, int lda, thrust::complex<T> const* B, int ldb,
            thrust::complex<T>* C, int ldc1, int ldc2)
{
    int i = blockIdx.x;
    int j = blockIdx.y*blockDim.x + threadIdx.x;
    int k = blockIdx.z;

    if( (i < ni) && (j<nj) && (k<nk) ) 
        C[ (k*ldc1 + i)*ldc2 + j ] = A[i*lda+j] * B[j*ldb+k];
}

template<typename T>
__global__ void kernel_element_wise_Aij_Bjk_Ckij( char transA, int ni, int nj, int nk,
            thrust::complex<T> const* A, int lda, thrust::complex<T> const* B, int ldb,
            thrust::complex<T>* C, int ldc1, int ldc2)
{
    int i = blockIdx.x;
    int j = blockIdx.y*blockDim.x + threadIdx.x;
    int k = blockIdx.z;

    if( (i < ni) && (j<nj) && (k<nk) ) {
        if(transA == 'N')
            C[ (k*ldc1 + i)*ldc2 + j ] = A[i*lda+j] * B[j*ldb+k];
        else if(transA == 'C')
            C[ (k*ldc1 + i)*ldc2 + j ] = conj(A[i*lda+j]) * B[j*ldb+k];
    }
}

// C[u][w] = alpha * sum_a A[u][w][a] * B[u][a]
void Auwn_Bun_Cuw(int nu, int nw, int na, std::complex<double> alpha, std::complex<double> const* A,
                        std::complex<double> const* B, std::complex<double> *C)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw)); 
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Auwn_Bun_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(A),
                                   reinterpret_cast<thrust::complex<double> const*>(B),
                                   reinterpret_cast<thrust::complex<double> *>(C));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


void Auwn_Bun_Cuw(int nu, int nw, int na, std::complex<float> alpha, std::complex<float> const* A,
                        std::complex<float> const* B, std::complex<float> *C)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw));
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Auwn_Bun_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(A),
                                   reinterpret_cast<thrust::complex<float> const*>(B),
                                   reinterpret_cast<thrust::complex<float> *>(C));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
void Awiu_Biu_Cuw(int nu, int nw, int na, std::complex<double> alpha, std::complex<double> const* A,
                        double const* B, int ldb, std::complex<double> *C, int ldc)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw));
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Awiu_Biu_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(A),
                                   B,ldb, 
                                   reinterpret_cast<thrust::complex<double> *>(C),ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


void Awiu_Biu_Cuw(int nu, int nw, int na, std::complex<float> alpha, std::complex<float> const* A,
                        float const* B, int ldb, std::complex<float> *C, int ldc)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw));
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Awiu_Biu_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(A),
                                   B,ldb,
                                   reinterpret_cast<thrust::complex<float> *>(C),ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void Awiu_Biu_Cuw(int nu, int nw, int na, std::complex<double> alpha, std::complex<double> const* A,
                        std::complex<double> const* B, int ldb, std::complex<double> *C, int ldc)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw));
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Awiu_Biu_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<double> const>(alpha),
                                   reinterpret_cast<thrust::complex<double> const*>(A),
                                   reinterpret_cast<thrust::complex<double> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<double> *>(C),ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


void Awiu_Biu_Cuw(int nu, int nw, int na, std::complex<float> alpha, std::complex<float> const* A,
                        std::complex<float> const* B, int ldb, std::complex<float> *C, int ldc)
{
  if(size_t(nw) > MAX_THREADS_PER_DIM) throw;
  size_t nthr = std::max(size_t(1),MAX_THREADS_PER_DIM/size_t(nw));
  size_t nbks = (nu + nthr - 1)/nthr;
  dim3 grid_dim(nbks,1,1);
  dim3 block_dim(nthr,nw,1);
  kernel_Awiu_Biu_Cuw<<<grid_dim,block_dim>>>(nu,nw,na,
                                   static_cast<thrust::complex<float> const>(alpha),
                                   reinterpret_cast<thrust::complex<float> const*>(A),
                                   reinterpret_cast<thrust::complex<float> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<float> *>(C),ldc);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

// element-wise C[k][i][j] = A[i][j] * B[j][k]
void element_wise_Aij_Bjk_Ckij(char transA, int ni, int nj, int nk,
            double const* A, int lda,
            std::complex<double> const* B, int ldb, std::complex<double> *C, int ldc1, int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1)/MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni,nbks,nk);
  kernel_element_wise_Aij_Bjk_Ckij<<<grid_dim,MAX_THREADS_PER_DIM>>>(transA,ni,nj,nk,
                                   A,lda, 
                                   reinterpret_cast<thrust::complex<double> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<double> *>(C),ldc1,ldc2);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA, int ni, int nj, int nk,
            float const* A, int lda,
            std::complex<float> const* B, int ldb, std::complex<float> *C, int ldc1, int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1)/MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni,nbks,nk);
  kernel_element_wise_Aij_Bjk_Ckij<<<grid_dim,MAX_THREADS_PER_DIM>>>(transA,ni,nj,nk,
                                   A,lda, 
                                   reinterpret_cast<thrust::complex<float> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<float> *>(C),ldc1,ldc2);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA, int ni, int nj, int nk, 
            std::complex<double> const* A, int lda,
            std::complex<double> const* B, int ldb, std::complex<double> *C, int ldc1, int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1)/MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni,nbks,nk);
  kernel_element_wise_Aij_Bjk_Ckij<<<grid_dim,MAX_THREADS_PER_DIM>>>(transA,ni,nj,nk,
                                   reinterpret_cast<thrust::complex<double> const*>(A),lda,
                                   reinterpret_cast<thrust::complex<double> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<double> *>(C),ldc1,ldc2);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA, int ni, int nj, int nk,
            std::complex<float> const* A, int lda,
            std::complex<float> const* B, int ldb, std::complex<float> *C, int ldc1, int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1)/MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni,nbks,nk);
  kernel_element_wise_Aij_Bjk_Ckij<<<grid_dim,MAX_THREADS_PER_DIM>>>(transA,ni,nj,nk,
                                   reinterpret_cast<thrust::complex<float> const*>(A),lda,
                                   reinterpret_cast<thrust::complex<float> const*>(B),ldb,
                                   reinterpret_cast<thrust::complex<float> *>(C),ldc1,ldc2);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


}
