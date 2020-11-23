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

#include <cassert>
#include <complex>
#include <cuda.h>
#include <thrust/complex.h>
#include <cuda_runtime.h>
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// very sloppy, needs improvement!!!!
// A[nocc_tot][nmo_tot][nwalk]
// B[Q][K][nocc_max][nmo_max][nwalk]
template<typename T, typename T2>
__global__ void kernel_KaKjw_to_QKajw(int nwalk,
                                      int nkpts,
                                      int npol,
                                      int nmo_max,
                                      int nmo_tot,
                                      int nocc_max,
                                      int* nmo,
                                      int* nmo0,
                                      int* nocc,
                                      int* nocc0,
                                      int* QKtok2,
                                      T const* A,
                                      T2* B)
{
  int Q   = blockIdx.x;
  int K   = blockIdx.y;
  int pol = blockIdx.z;
  if (Q >= nkpts || K >= nkpts || pol > npol)
    return;
  int QK  = QKtok2[Q * nkpts + K];
  int na0 = nocc0[K];
  int nj0 = nmo0[QK];
  int na  = nocc[K];
  int nj  = nmo[QK];

  T const* A_(A + (na0 * npol * nmo_tot + nj0) * nwalk);
  T2* B_(B + ((Q * nkpts + K) * nocc_max) * npol * nmo_max * nwalk);

  if (threadIdx.x >= nj)
    return;
  if (threadIdx.y >= nwalk)
    return;

  for (int a = 0, a0 = pol * nmo_max * nwalk, a1 = pol * nmo_tot * nwalk; a < na;
       a++, a0 += npol * nmo_max * nwalk, a1 += npol * nmo_tot * nwalk)
    for (int j = threadIdx.x; j < nj; j += blockDim.x)
      for (int n = threadIdx.y; n < nwalk; n += blockDim.y)
        B_[a0 + j * nwalk + n] = static_cast<T2>(A_[a1 + j * nwalk + n]);
}

template<typename T, typename T2>
__global__ void kernel_KaKjw_to_QKajw(int nwalk,
                                      int nkpts,
                                      int npol,
                                      int nmo_max,
                                      int nmo_tot,
                                      int nocc_max,
                                      int* nmo,
                                      int* nmo0,
                                      int* nocc,
                                      int* nocc0,
                                      int* QKtok2,
                                      thrust::complex<T> const* A,
                                      thrust::complex<T2>* B)
{
  int Q   = blockIdx.x;
  int K   = blockIdx.y;
  int pol = blockIdx.z;
  if (Q >= nkpts || K >= nkpts || pol > npol)
    return;
  int QK  = QKtok2[Q * nkpts + K];
  int na0 = nocc0[K];
  int nj0 = nmo0[QK];
  int na  = nocc[K];
  int nj  = nmo[QK];

  thrust::complex<T> const* A_(A + (na0 * npol * nmo_tot + nj0) * nwalk);
  thrust::complex<T2>* B_(B + ((Q * nkpts + K) * nocc_max) * npol * nmo_max * nwalk);

  if (threadIdx.x >= nj)
    return;
  if (threadIdx.y >= nwalk)
    return;

  for (int a = 0, a0 = pol * nmo_max * nwalk, a1 = pol * nmo_tot * nwalk; a < na;
       a++, a0 += npol * nmo_max * nwalk, a1 += npol * nmo_tot * nwalk)
  {
    for (int j = threadIdx.x; j < nj; j += blockDim.x)
      for (int n = threadIdx.y; n < nwalk; n += blockDim.y)
        B_[a0 + j * nwalk + n] = static_cast<thrust::complex<T2>>(A_[a1 + j * nwalk + n]);
  }
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    double const* A,
                    double* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, A, B);
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    float const* A,
                    float* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, A, B);
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    double const* A,
                    float* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, A, B);
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    std::complex<float> const* A,
                    std::complex<float>* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, reinterpret_cast<thrust::complex<float> const*>(A),
                                                 reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    std::complex<double> const* A,
                    std::complex<double>* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, reinterpret_cast<thrust::complex<double> const*>(A),
                                                 reinterpret_cast<thrust::complex<double>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

void KaKjw_to_QKajw(int nwalk,
                    int nkpts,
                    int npol,
                    int nmo_max,
                    int nmo_tot,
                    int nocc_max,
                    int* nmo,
                    int* nmo0,
                    int* nocc,
                    int* nocc0,
                    int* QKtok2,
                    std::complex<double> const* A,
                    std::complex<float>* B)
{
  int xblock_dim = 16;
  int yblock_dim = std::min(nwalk, 32);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, nkpts, npol);
  kernel_KaKjw_to_QKajw<<<grid_dim, block_dim>>>(nwalk, nkpts, npol, nmo_max, nmo_tot, nocc_max, nmo, nmo0, nocc, nocc0,
                                                 QKtok2, reinterpret_cast<thrust::complex<double> const*>(A),
                                                 reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError(), "KaKjw_to_QKajw");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "KaKjw_to_QKajw");
}

} // namespace kernels
