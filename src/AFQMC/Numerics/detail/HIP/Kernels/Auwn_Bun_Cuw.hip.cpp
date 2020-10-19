///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <complex>
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include "uninitialized_array.hpp"
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
// C[u][w] = alpha * sum_a A[u][w][a] * B[u][a]
template<typename T>
__global__ void kernel_Auwn_Bun_Cuw(int nu,
                                    int nw,
                                    int na,
                                    thrust::complex<T> const alpha,
                                    thrust::complex<T> const* A,
                                    thrust::complex<T> const* B,
                                    thrust::complex<T>* C)
{
  int nu_per_block = blockDim.x;
  int u            = blockIdx.x * nu_per_block + threadIdx.x;
  int w            = threadIdx.y;

  if ((u < nu) && (w < nw))
  {
    thrust::complex<T> Cuw = 0;
    thrust::complex<T> const* A_(A + (u * nw + w) * na);
    thrust::complex<T> const* B_(B + u * na);
    for (int a = 0; a < na; ++a, ++A_, ++B_)
      Cuw += (*A_) * (*B_);
    C[u * nw + w] = alpha * Cuw;
  }
}

// NOT OPTIMAL: poor memory access pattern
template<typename T>
__global__ void kernel_Awiu_Biu_Cuw(int nu,
                                    int nw,
                                    int ni,
                                    thrust::complex<T> const alpha,
                                    thrust::complex<T> const* A,
                                    T const* B,
                                    int ldb,
                                    thrust::complex<T>* C,
                                    int ldc)
{
  int nu_per_block = blockDim.x;
  int u            = blockIdx.x * nu_per_block + threadIdx.x;
  int w            = threadIdx.y;

  if ((u < nu) && (w < nw))
  {
    thrust::complex<T> Cuw = 0;
    thrust::complex<T> const* A_(A + w * ni * nu + u);
    T const* B_(B + u);
    for (int i = 0; i < ni; ++i, A_ += nu, B_ += ldb)
      Cuw += (*A_) * (*B_);
    C[u * ldc + w] = alpha * Cuw;
  }
}

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
template<typename T>
__global__ void kernel_Awiu_Biu_Cuw(int nu,
                                    int nw,
                                    int ni,
                                    thrust::complex<T> const alpha,
                                    thrust::complex<T> const* A,
                                    thrust::complex<T> const* B,
                                    int ldb,
                                    thrust::complex<T>* C,
                                    int ldc)
{
  int nu_per_block = blockDim.x;
  int u            = blockIdx.x * nu_per_block + threadIdx.x;
  int w            = threadIdx.y;

  if ((u < nu) && (w < nw))
  {
    thrust::complex<T> Cuw = 0;
    thrust::complex<T> const* A_(A + w * ni * nu + u);
    thrust::complex<T> const* B_(B + u);
    for (int i = 0; i < ni; ++i, A_ += nu, B_ += ldb)
      Cuw += (*A_) * (*B_);
    C[u * ldc + w] = alpha * Cuw;
  }
}

// Cik = sum_j Aijk * Bkj  ( nthreads per block 32 )
template<typename T, typename T1>
__global__ void kernel_Aijk_Bkj_Cik(int ni,
                                    int nj,
                                    int nk,
                                    thrust::complex<T> const* A,
                                    int lda,
                                    int stride,
                                    T1 const* B,
                                    int ldb,
                                    thrust::complex<T>* C,
                                    int ldc)
{
  __shared__ uninitialized_array<thrust::complex<T>, 32> cache;
  int k = blockIdx.x;
  int i = blockIdx.y;
  if ((i < ni) && (k < nk))
  {
    cache[threadIdx.x] = thrust::complex<T>(0.0);
    auto A_(A + i * stride + k + threadIdx.x * lda);
    auto B_(B + k * ldb + threadIdx.x);
    auto Bn_(B + k * ldb + nj);
    while (B_ < Bn_)
    {
      cache[threadIdx.x] += (*A_) * static_cast<thrust::complex<T>>(*B_);
      A_ += blockDim.x * lda;
      B_ += blockDim.x;
    }

    __syncthreads();
    int j = 16;
    while (j > 0)
    {
      if (threadIdx.x < j)
        cache[threadIdx.x] += cache[threadIdx.x + j];
      __syncthreads();
      j /= 2; //not sure bitwise operations are actually faster
    }
    if (threadIdx.x == 0)
      *(C + i * ldc + k) += cache[0];
  }
}

// A[w][i][j] = B[i][w][j]
template<typename T, typename T1>
__global__ void kernel_viwj_vwij(int nw, int ni, int i0, int iN, thrust::complex<T> const* B, thrust::complex<T1>* A)
{
  int w = blockIdx.x;
  int i = blockIdx.y + i0;
  if ((w < nw) && (i < ni))
  {
    int j = threadIdx.x;
    auto A_(A + (w * ni + i) * ni);
    auto B_(B + (i * nw + w) * ni);
    while (j < ni)
    {
      A_[j] = static_cast<thrust::complex<T1>>(B_[j]);
      j += blockDim.x;
    }
  }
}

// element-wise C[k][i][j] = A[i][j] * B[j][k]
template<typename T>
__global__ void kernel_element_wise_Aij_Bjk_Ckij(char transA,
                                                 int ni,
                                                 int nj,
                                                 int nk,
                                                 T const* A,
                                                 int lda,
                                                 thrust::complex<T> const* B,
                                                 int ldb,
                                                 thrust::complex<T>* C,
                                                 int ldc1,
                                                 int ldc2)
{
  int i = blockIdx.x;
  int j = blockIdx.y * blockDim.x + threadIdx.x;
  int k = blockIdx.z;

  if ((i < ni) && (j < nj) && (k < nk))
    C[(k * ldc1 + i) * ldc2 + j] = A[i * lda + j] * B[j * ldb + k];
}

template<typename T>
__global__ void kernel_element_wise_Aij_Bjk_Ckij(char transA,
                                                 int ni,
                                                 int nj,
                                                 int nk,
                                                 thrust::complex<T> const* A,
                                                 int lda,
                                                 thrust::complex<T> const* B,
                                                 int ldb,
                                                 thrust::complex<T>* C,
                                                 int ldc1,
                                                 int ldc2)
{
  int i = blockIdx.x;
  int j = blockIdx.y * blockDim.x + threadIdx.x;
  int k = blockIdx.z;

  if ((i < ni) && (j < nj) && (k < nk))
  {
    if (transA == 'N')
      C[(k * ldc1 + i) * ldc2 + j] = A[i * lda + j] * B[j * ldb + k];
    else if (transA == 'C')
      C[(k * ldc1 + i) * ldc2 + j] = conj(A[i * lda + j]) * B[j * ldb + k];
  }
}

// Ckji = Aij * Bjk
template<typename T, typename T2>
__global__ void kernel_element_wise_Aij_Bjk_Ckji(int ni,
                                                 int nj,
                                                 int nk,
                                                 T2 const* A,
                                                 int lda,
                                                 thrust::complex<T> const* B,
                                                 int ldb,
                                                 thrust::complex<T>* C,
                                                 int ldc,
                                                 int stride)
{
  // hard-coded to TILE_DIM=32
  int TILE_DIM = 32;
  __shared__ uninitialized_array<T2, 32 * 32> Acache;
  __shared__ uninitialized_array<thrust::complex<T>, 32> Bcache;

  int k = blockIdx.z;
  int j = blockIdx.x * TILE_DIM + threadIdx.x;
  int i = blockIdx.y * TILE_DIM + threadIdx.y;

  if ((k < nk) && (j < nj))
  {
    int n(threadIdx.y);
    while ((i < ni) && (n < TILE_DIM))
    {
      Acache[n * 32 + threadIdx.x] = A[i * lda + j];
      n += blockDim.y;
      i += blockDim.y;
    }
    if (threadIdx.y == 0)
      Bcache[threadIdx.x] = B[j * ldb + k];
  }

  __syncthreads();

  // subtle interchange of threadIdx.x/threadIdx.y
  i = blockIdx.y * TILE_DIM + threadIdx.x;
  j = blockIdx.x * TILE_DIM + threadIdx.y;

  if ((k < nk) && (i < ni))
  {
    int n(threadIdx.y);
    while ((j < nj) && (n < TILE_DIM))
    {
      C[k * stride + j * ldc + i] = static_cast<thrust::complex<T>>(Acache[threadIdx.x * 32 + n]) * Bcache[n];
      n += blockDim.y;
      j += blockDim.y;
    }
  }
}

// C[u][w] = alpha * sum_a A[u][w][a] * B[u][a]
void Auwn_Bun_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  std::complex<double> const* B,
                  std::complex<double>* C)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Auwn_Bun_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A),
                     reinterpret_cast<thrust::complex<double> const*>(B),
                     reinterpret_cast<thrust::complex<double>*>(C));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


void Auwn_Bun_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  std::complex<float> const* B,
                  std::complex<float>* C)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Auwn_Bun_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A),
                     reinterpret_cast<thrust::complex<float> const*>(B), reinterpret_cast<thrust::complex<float>*>(C));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// C[u][w] = alpha * sum_i A[w][i][u] * B[i][u]
void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  double const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Awiu_Biu_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A), B, ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  float const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Awiu_Biu_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A), B, ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<double> alpha,
                  std::complex<double> const* A,
                  std::complex<double> const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Awiu_Biu_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(A),
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


void Awiu_Biu_Cuw(int nu,
                  int nw,
                  int na,
                  std::complex<float> alpha,
                  std::complex<float> const* A,
                  std::complex<float> const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc)
{
  if (size_t(nw) > MAX_THREADS_PER_DIM)
    throw;
  size_t nthr = std::max(size_t(1), MAX_THREADS_PER_DIM / size_t(nw));
  size_t nbks = (nu + nthr - 1) / nthr;
  dim3 grid_dim(nbks, 1, 1);
  dim3 block_dim(nthr, nw, 1);
  hipLaunchKernelGGL(kernel_Awiu_Biu_Cuw, dim3(grid_dim), dim3(block_dim), 0, 0, nu, nw, na,
                     static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(A),
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


// C[i][k] = sum_i A[i][j][k] * B[k][j]
void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<double> const* A,
                  int lda,
                  int stride,
                  std::complex<double> const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc)
{
  // expect nk >> ni,nj
  dim3 grid_dim(nk, ni, 1);
  hipLaunchKernelGGL(kernel_Aijk_Bkj_Cik, dim3(grid_dim), dim3(32), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<double> const*>(A), lda, stride,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<double> const* A,
                  int lda,
                  int stride,
                  double const* B,
                  int ldb,
                  std::complex<double>* C,
                  int ldc)
{
  // expect nk >> ni,nj
  dim3 grid_dim(nk, ni, 1);
  hipLaunchKernelGGL(kernel_Aijk_Bkj_Cik, dim3(grid_dim), dim3(32), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<double> const*>(A), lda, stride, B, ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<float> const* A,
                  int lda,
                  int stride,
                  std::complex<float> const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc)
{
  // expect nk >> ni,nj
  dim3 grid_dim(nk, ni, 1);
  hipLaunchKernelGGL(kernel_Aijk_Bkj_Cik, dim3(grid_dim), dim3(32), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<float> const*>(A), lda, stride,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void Aijk_Bkj_Cik(int ni,
                  int nj,
                  int nk,
                  std::complex<float> const* A,
                  int lda,
                  int stride,
                  float const* B,
                  int ldb,
                  std::complex<float>* C,
                  int ldc)
{
  // expect nk >> ni,nj
  dim3 grid_dim(nk, ni, 1);
  hipLaunchKernelGGL(kernel_Aijk_Bkj_Cik, dim3(grid_dim), dim3(32), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<float> const*>(A), lda, stride, B, ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// v[w][i][j] = v[i][w][j]
void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<double> const* B, std::complex<double>* A)
{
  // expect ni > nw
  dim3 grid_dim(nw, (iN - i0), 1);
  hipLaunchKernelGGL(kernel_viwj_vwij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, nw, ni, i0, iN,
                     reinterpret_cast<thrust::complex<double> const*>(B),
                     reinterpret_cast<thrust::complex<double>*>(A));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<double> const* B, std::complex<float>* A)
{
  // expect ni > nw
  dim3 grid_dim(nw, (iN - i0), 1);
  hipLaunchKernelGGL(kernel_viwj_vwij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, nw, ni, i0, iN,
                     reinterpret_cast<thrust::complex<double> const*>(B), reinterpret_cast<thrust::complex<float>*>(A));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<float> const* B, std::complex<double>* A)
{
  // expect ni > nw
  dim3 grid_dim(nw, (iN - i0), 1);
  hipLaunchKernelGGL(kernel_viwj_vwij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, nw, ni, i0, iN,
                     reinterpret_cast<thrust::complex<float> const*>(B), reinterpret_cast<thrust::complex<double>*>(A));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void viwj_vwij(int nw, int ni, int i0, int iN, std::complex<float> const* B, std::complex<float>* A)
{
  // expect ni > nw
  dim3 grid_dim(nw, (iN - i0), 1);
  hipLaunchKernelGGL(kernel_viwj_vwij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, nw, ni, i0, iN,
                     reinterpret_cast<thrust::complex<float> const*>(B), reinterpret_cast<thrust::complex<float>*>(A));
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

// element-wise C[k][i][j] = A[i][j] * B[j][k]
void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               double const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc1,
                               int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni, nbks, nk);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, transA, ni, nj,
                     nk, A, lda, reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc1, ldc2);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               float const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc1,
                               int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni, nbks, nk);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, transA, ni, nj,
                     nk, A, lda, reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc1, ldc2);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               std::complex<double> const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc1,
                               int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni, nbks, nk);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, transA, ni, nj,
                     nk, reinterpret_cast<thrust::complex<double> const*>(A), lda,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc1, ldc2);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckij(char transA,
                               int ni,
                               int nj,
                               int nk,
                               std::complex<float> const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc1,
                               int ldc2)
{
  // setup for nj >> ni,nk
  size_t nbks = (nj + MAX_THREADS_PER_DIM - 1) / MAX_THREADS_PER_DIM;
  dim3 grid_dim(ni, nbks, nk);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckij, dim3(grid_dim), dim3(MAX_THREADS_PER_DIM), 0, 0, transA, ni, nj,
                     nk, reinterpret_cast<thrust::complex<float> const*>(A), lda,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc1, ldc2);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


// element-wise C[k][j][i] = A[i][j] * B[j][k]
void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               double const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc,
                               int stride)
{
  // setup for nj >> ni,nk
  size_t nthr  = 32;
  size_t nthrj = 8;
  size_t ib    = (ni + nthr - 1) / nthr;
  size_t jb    = (nj + nthr - 1) / nthr;
  dim3 grid_dim(jb, ib, nk);
  dim3 block_dim(nthr, nthrj, 1);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckji, dim3(grid_dim), dim3(block_dim), 0, 0, ni, nj, nk, A, lda,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc, stride);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               std::complex<double> const* A,
                               int lda,
                               std::complex<double> const* B,
                               int ldb,
                               std::complex<double>* C,
                               int ldc,
                               int stride)
{
  // setup for nj >> ni,nk
  size_t nthr  = 32;
  size_t nthrj = 8;
  size_t ib    = (ni + nthr - 1) / nthr;
  size_t jb    = (nj + nthr - 1) / nthr;
  // jb goes along x since this is the fast index in Aij, needed for better memory access patterns
  dim3 grid_dim(jb, ib, nk);
  dim3 block_dim(nthr, nthrj, 1);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckji, dim3(grid_dim), dim3(block_dim), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<double> const*>(A), lda,
                     reinterpret_cast<thrust::complex<double> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<double>*>(C), ldc, stride);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               float const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc,
                               int stride)
{
  // setup for nj >> ni,nk
  size_t nthr  = 32;
  size_t nthrj = 8;
  size_t ib    = (ni + nthr - 1) / nthr;
  size_t jb    = (nj + nthr - 1) / nthr;
  dim3 grid_dim(jb, ib, nk);
  dim3 block_dim(nthr, nthrj, 1);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckji, dim3(grid_dim), dim3(block_dim), 0, 0, ni, nj, nk, A, lda,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc, stride);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}

void element_wise_Aij_Bjk_Ckji(int ni,
                               int nj,
                               int nk,
                               std::complex<float> const* A,
                               int lda,
                               std::complex<float> const* B,
                               int ldb,
                               std::complex<float>* C,
                               int ldc,
                               int stride)
{
  // setup for nj >> ni,nk
  size_t nthr  = 32;
  size_t nthrj = 8;
  size_t ib    = (ni + nthr - 1) / nthr;
  size_t jb    = (nj + nthr - 1) / nthr;
  dim3 grid_dim(jb, ib, nk);
  dim3 block_dim(nthr, nthrj, 1);
  hipLaunchKernelGGL(kernel_element_wise_Aij_Bjk_Ckji, dim3(grid_dim), dim3(block_dim), 0, 0, ni, nj, nk,
                     reinterpret_cast<thrust::complex<float> const*>(A), lda,
                     reinterpret_cast<thrust::complex<float> const*>(B), ldb,
                     reinterpret_cast<thrust::complex<float>*>(C), ldc, stride);
  qmc_hip::hip_kernel_check(hipGetLastError());
  qmc_hip::hip_kernel_check(hipDeviceSynchronize());
}


} // namespace kernels
