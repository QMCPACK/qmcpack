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
template<typename T, typename T2>
__global__ void kernel_ajw_to_waj(int na, int nj, int nw, int inca, T const* A, T2* B)
{
  int a = blockIdx.x;
  if (a >= na)
    return;
  T const* A_(A + inca * a);
  int i   = threadIdx.x;
  int njw = nj * nw;
  while (i < njw)
  {
    int j                    = i / nw;
    int w                    = i % nw;
    B[(w * na + a) * nj + j] = static_cast<T2>(A_[i]);
    i += blockDim.x;
  }
}

template<typename T, typename T2>
__global__ void kernel_ajw_to_waj(int na, int nj, int nw, int inca, thrust::complex<T> const* A, thrust::complex<T2>* B)
{
  int a = blockIdx.x;
  if (a >= na)
    return;
  thrust::complex<T> const* A_(A + inca * a);
  int i   = threadIdx.x;
  int njw = nj * nw;
  while (i < njw)
  {
    int j                    = i / nw;
    int w                    = i % nw;
    B[(w * na + a) * nj + j] = static_cast<thrust::complex<T2>>(A_[i]);
    i += blockDim.x;
  }
}

template<typename T, typename T2>
__global__ void kernel_transpose_wabn_to_wban(int nwalk,
                                              int na,
                                              int nb,
                                              int nchol,
                                              thrust::complex<T> const* Tab,
                                              thrust::complex<T2>* Tba)
{
  int w = blockIdx.x;
  int a = blockIdx.y;
  int b = blockIdx.z;

  thrust::complex<T> const* A_(Tab + ((w * na + a) * nb + b) * nchol);
  thrust::complex<T2>* B_(Tba + ((w * nb + b) * na + a) * nchol);
  int i = threadIdx.x;
  while (i < nchol)
  {
    B_[i] = static_cast<thrust::complex<T2>>(A_[i]);
    i += blockDim.x;
  }
}

void ajw_to_waj(int na, int nj, int nw, int inca, double const* A, double* B)
{
  kernel_ajw_to_waj<<<na, 128>>>(na, nj, nw, inca, A, B);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, float const* A, float* B)
{
  kernel_ajw_to_waj<<<na, 128>>>(na, nj, nw, inca, A, B);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<double> const* A, std::complex<double>* B)
{
  kernel_ajw_to_waj<<<na, 128>>>(na, nj, nw, inca, reinterpret_cast<thrust::complex<double> const*>(A),
                                 reinterpret_cast<thrust::complex<double>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void ajw_to_waj(int na, int nj, int nw, int inca, std::complex<float> const* A, std::complex<float>* B)
{
  kernel_ajw_to_waj<<<na, 128>>>(na, nj, nw, inca, reinterpret_cast<thrust::complex<float> const*>(A),
                                 reinterpret_cast<thrust::complex<float>*>(B));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<double>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  kernel_transpose_wabn_to_wban<<<grid_dim, 32>>>(nwalk, na, nb, nchol,
                                                  reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                  reinterpret_cast<thrust::complex<double>*>(Tba));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<float>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  kernel_transpose_wabn_to_wban<<<grid_dim, 32>>>(nwalk, na, nb, nchol,
                                                  reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                  reinterpret_cast<thrust::complex<float>*>(Tba));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<double> const* Tab,
                            std::complex<float>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  kernel_transpose_wabn_to_wban<<<grid_dim, 32>>>(nwalk, na, nb, nchol,
                                                  reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                  reinterpret_cast<thrust::complex<float>*>(Tba));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void transpose_wabn_to_wban(int nwalk,
                            int na,
                            int nb,
                            int nchol,
                            std::complex<float> const* Tab,
                            std::complex<double>* Tba)
{
  dim3 grid_dim(nwalk, na, nb);
  kernel_transpose_wabn_to_wban<<<grid_dim, 32>>>(nwalk, na, nb, nchol,
                                                  reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                  reinterpret_cast<thrust::complex<double>*>(Tba));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
