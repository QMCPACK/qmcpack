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
#include <thrust/system/cuda/detail/core/util.h>
#include "AFQMC/Numerics/detail/CUDA/Kernels/cuda_settings.h"
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// X[m,ni,iw] = rand[m,ni,iw] + im * ( vbias[m,iw] - vMF[m]  )
// HW[ni,iw] = sum_m [ im * ( vMF[m] - vbias[m,iw] ) *
//                     ( rand[m,ni,iw] + halfim * ( vbias[m,iw] - vMF[m] ) ) ]
//           = sum_m [ im * ( vMF[m] - vbias[m,iw] ) *
//                     ( X[m,ni,iw] - halfim * ( vbias[m,iw] - vMF[m] ) ) ]
// MF[ni,iw] = sum_m ( im * X[m,ni,iw] * vMF[m] )
template<typename T, typename T2, typename T3>
__global__ void kernel_construct_X(int nCV,
                                   int nsteps,
                                   int nwalk,
                                   T sqrtdt,
                                   T vbound,
                                   thrust::complex<T> const* vMF,
                                   thrust::complex<T2> const* vbias,
                                   thrust::complex<T>* HW,
                                   thrust::complex<T>* MF,
                                   thrust::complex<T3>* X)
{
  int ni = blockIdx.x;
  int iw = blockIdx.y;
  if (ni >= nsteps || iw >= nwalk)
    return;
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 2 * REDUCE_BLOCK_SIZE> cache;
  cache[2 * threadIdx.x]     = thrust::complex<T>(0.0);
  cache[2 * threadIdx.x + 1] = thrust::complex<T>(0.0);

  thrust::complex<T> im(0.0, 1.0);
  thrust::complex<T> const* vmf_(vMF);
  thrust::complex<T2> const* vb_(vbias + iw);
  thrust::complex<T3>* X_(X + ni * nwalk + iw);
  int dvb_ = nwalk;
  int dX   = nsteps * nwalk;

  thrust::complex<T> vmf_t;
  thrust::complex<T> vb_t;
  for (int m = threadIdx.x; m < nCV; m += blockDim.x)
  {
    if (abs(vmf_[m]) > vbound)
      vmf_t = sqrtdt * vmf_[m] / (abs(vmf_[m]) / vbound);
    else
      vmf_t = sqrtdt * vmf_[m];
    if (abs(vb_[m * dvb_]) > static_cast<T2>(vbound * sqrtdt))
      vb_t = static_cast<thrust::complex<T>>(vb_[m * dvb_] / (abs(vb_[m * dvb_]) / static_cast<T2>(vbound * sqrtdt)));
    else
      vb_t = static_cast<thrust::complex<T>>(vb_[m * dvb_]);

    thrust::complex<T> vdiff = (im * (vb_t - vmf_t));
    X_[m * dX] += static_cast<thrust::complex<T3>>(vdiff);
    cache[2 * threadIdx.x] -= vdiff * (static_cast<thrust::complex<T>>(X_[m * dX]) - 0.5 * vdiff);
    cache[2 * threadIdx.x + 1] += im * static_cast<thrust::complex<T>>(X_[m * dX]) * (sqrtdt * vmf_[m]);
  }

  __syncthreads(); // required because later on the current thread is accessing
                   // data written by another thread
  int i = REDUCE_BLOCK_SIZE / 2;
  while (i > 0)
  {
    if (threadIdx.x < i)
    {
      cache[2 * threadIdx.x] += cache[2 * (threadIdx.x + i)];
      cache[2 * threadIdx.x + 1] += cache[2 * (threadIdx.x + i) + 1];
    }
    __syncthreads();
    i /= 2; //not sure bitwise operations are actually faster
  }
  if (threadIdx.x == 0)
  {
    HW[ni * nwalk + iw] = cache[0];
    MF[ni * nwalk + iw] = cache[1];
  }
}

template<typename T, typename T1>
__global__ void kernel_construct_X_free_projection(int nCV,
                                                   int nsteps,
                                                   int nwalk,
                                                   T sqrtdt,
                                                   thrust::complex<T> const* vMF,
                                                   thrust::complex<T>* MF,
                                                   thrust::complex<T1>* X)
{
  int ni = blockIdx.x;
  int iw = blockIdx.y;
  if (ni >= nsteps || iw >= nwalk)
    return;
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, REDUCE_BLOCK_SIZE> cache;
  cache[threadIdx.x] = thrust::complex<T>(0.0);

  thrust::complex<T> im(0.0, 1.0);
  thrust::complex<T> const* vmf_(vMF);
  thrust::complex<T1>* X_(X + ni * nwalk + iw);
  int dX = nsteps * nwalk;

  for (int m = threadIdx.x; m < nCV; m += blockDim.x)
    cache[threadIdx.x] += im * static_cast<thrust::complex<T>>(X_[m * dX]) * vmf_[m] * sqrtdt;

  __syncthreads(); // required because later on the current thread is accessing
                   // data written by another thread
  int i = REDUCE_BLOCK_SIZE / 2;
  while (i > 0)
  {
    if (threadIdx.x < i)
      cache[threadIdx.x] += cache[threadIdx.x + i];
    __syncthreads();
    i /= 2; //not sure bitwise operations are actually faster
  }
  if (threadIdx.x == 0)
  {
    MF[ni * nwalk + iw] = cache[0];
  }
}

void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<double> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<double>* X)
{
  dim3 block_dim(REDUCE_BLOCK_SIZE, 1, 1);
  dim3 grid_dim(nsteps, nwalk, 1);
  if (free_projection)
    kernel_construct_X_free_projection<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt,
                                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                                reinterpret_cast<thrust::complex<double>*>(X));
  else
    kernel_construct_X<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt, vbound,
                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                reinterpret_cast<thrust::complex<double> const*>(vbias),
                                                reinterpret_cast<thrust::complex<double>*>(HW),
                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                reinterpret_cast<thrust::complex<double>*>(X));
  qmc_cuda::cuda_check(cudaGetLastError(), "construct_X");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "construct_X");
}
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<float> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<float>* X)
{
  dim3 block_dim(REDUCE_BLOCK_SIZE, 1, 1);
  dim3 grid_dim(nsteps, nwalk, 1);
  if (free_projection)
    kernel_construct_X_free_projection<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt,
                                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                                reinterpret_cast<thrust::complex<float>*>(X));
  else
    kernel_construct_X<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt, vbound,
                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                reinterpret_cast<thrust::complex<float> const*>(vbias),
                                                reinterpret_cast<thrust::complex<double>*>(HW),
                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                reinterpret_cast<thrust::complex<float>*>(X));
  qmc_cuda::cuda_check(cudaGetLastError(), "construct_X");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "construct_X");
}
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<double> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<float>* X)
{
  dim3 block_dim(REDUCE_BLOCK_SIZE, 1, 1);
  dim3 grid_dim(nsteps, nwalk, 1);
  if (free_projection)
    kernel_construct_X_free_projection<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt,
                                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                                reinterpret_cast<thrust::complex<float>*>(X));
  else
    kernel_construct_X<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt, vbound,
                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                reinterpret_cast<thrust::complex<double> const*>(vbias),
                                                reinterpret_cast<thrust::complex<double>*>(HW),
                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                reinterpret_cast<thrust::complex<float>*>(X));
  qmc_cuda::cuda_check(cudaGetLastError(), "construct_X");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "construct_X");
}
void construct_X(int nCV,
                 int nsteps,
                 int nwalk,
                 bool free_projection,
                 double sqrtdt,
                 double vbound,
                 std::complex<double> const* vMF,
                 std::complex<float> const* vbias,
                 std::complex<double>* HW,
                 std::complex<double>* MF,
                 std::complex<double>* X)
{
  dim3 block_dim(REDUCE_BLOCK_SIZE, 1, 1);
  dim3 grid_dim(nsteps, nwalk, 1);
  if (free_projection)
    kernel_construct_X_free_projection<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt,
                                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                                reinterpret_cast<thrust::complex<double>*>(X));
  else
    kernel_construct_X<<<grid_dim, block_dim>>>(nCV, nsteps, nwalk, sqrtdt, vbound,
                                                reinterpret_cast<thrust::complex<double> const*>(vMF),
                                                reinterpret_cast<thrust::complex<float> const*>(vbias),
                                                reinterpret_cast<thrust::complex<double>*>(HW),
                                                reinterpret_cast<thrust::complex<double>*>(MF),
                                                reinterpret_cast<thrust::complex<double>*>(X));
  qmc_cuda::cuda_check(cudaGetLastError(), "construct_X");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "construct_X");
}

} // namespace kernels
