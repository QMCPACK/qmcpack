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
#include <cuda_runtime.h>
#include "math.h"
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <cuda_runtime.h>
#include <thrust/system/cuda/detail/core/util.h>
#define ENABLE_CUDA 1
#include "AFQMC/Memory/CUDA/cuda_utilities.h"

namespace kernels
{
// Meant to be run with 1 block
template<class T>
__global__ void kernel_determinant_from_getrf(int N, T const* m, int lda, int const* piv, T LogOverlapFactor, T* det)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 256> tmp;
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 256> sg;
  int t = threadIdx.x;

  tmp[t] = T(0.0);
  sg[t]  = T(1.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m[ip * lda + ip] < 0.0)
    {
      tmp[t] += log(-m[ip * lda + ip]);
      sg[t] *= T(-1.0);
    }
    else
      tmp[t] += log(m[ip * lda + ip]);
    if (piv[ip] != (ip + 1))
      sg[t] *= T(-1.0);
  }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
    {
      tmp[0] += tmp[i];
      sg[0] *= sg[i];
    }
    *det = sg[0] * exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<class T>
__global__ void kernel_determinant_from_getrf(int N,
                                              thrust::complex<T> const* m,
                                              int lda,
                                              int const* piv,
                                              thrust::complex<T> LogOverlapFactor,
                                              thrust::complex<T>* det)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 256> tmp;
  int t = threadIdx.x;

  tmp[t] = thrust::complex<T>(0.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
    if (piv[ip] == (ip + 1))
    {
      tmp[t] += log(m[ip * lda + ip]);
    }
    else
    {
      tmp[t] += log(-m[ip * lda + ip]);
    }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
      tmp[0] += tmp[i];
    *det = exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<class T>
__global__ void kernel_strided_determinant_from_getrf(int N,
                                                      T const* m,
                                                      int lda,
                                                      int mstride,
                                                      int const* piv,
                                                      int pstride,
                                                      T LogOverlapFactor,
                                                      T* det,
                                                      int nbatch)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 64> tmp;
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 64> sg;
  int t     = threadIdx.x;
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;

  int const* piv_ = piv + batch * pstride;
  T const* m_     = m + batch * mstride;

  tmp[t] = T(0.0);
  sg[t]  = T(1.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m_[ip * lda + ip] < 0.0)
    {
      tmp[t] += log(-m_[ip * lda + ip]);
      sg[t] *= T(-1.0);
    }
    else
      tmp[t] += log(m_[ip * lda + ip]);
    if (piv_[ip] != (ip + 1))
      sg[t] *= T(-1.0);
  }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
    {
      tmp[0] += tmp[i];
      sg[0] *= sg[i];
    }
    *(det + batch) = sg[0] * exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<class T>
__global__ void kernel_strided_determinant_from_getrf(int N,
                                                      thrust::complex<T> const* m,
                                                      int lda,
                                                      int mstride,
                                                      int const* piv,
                                                      int pstride,
                                                      thrust::complex<T> LogOverlapFactor,
                                                      thrust::complex<T>* det,
                                                      int nbatch)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 64> tmp;
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;
  int t = threadIdx.x;

  tmp[t] = thrust::complex<T>(0.0);

  int const* piv_              = piv + batch * pstride;
  thrust::complex<T> const* m_ = m + batch * mstride;

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
    if (piv_[ip] == (ip + 1))
    {
      tmp[t] += log(m_[ip * lda + ip]);
    }
    else
    {
      tmp[t] += log(-m_[ip * lda + ip]);
    }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
      tmp[0] += tmp[i];
    *(det + batch) = exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<class T>
__global__ void kernel_batched_determinant_from_getrf(int N,
                                                      T* const* m,
                                                      int lda,
                                                      int const* piv,
                                                      int pstride,
                                                      T LogOverlapFactor,
                                                      T* det,
                                                      int nbatch)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 64> tmp;
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 64> sg;
  int t     = threadIdx.x;
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;

  int const* piv_ = piv + batch * pstride;
  T const* m_     = m[batch];

  tmp[t] = T(0.0);
  sg[t]  = T(1.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m_[ip * lda + ip] < 0.0)
    {
      tmp[t] += log(-m_[ip * lda + ip]);
      sg[t] *= T(-1.0);
    }
    else
      tmp[t] += log(m_[ip * lda + ip]);
    if (piv_[ip] != (ip + 1))
      sg[t] *= T(-1.0);
  }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
    {
      tmp[0] += tmp[i];
      sg[0] *= sg[i];
    }
    *(det + batch) = sg[0] * exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<class T>
__global__ void kernel_batched_determinant_from_getrf(int N,
                                                      thrust::complex<T>* const* m,
                                                      int lda,
                                                      int const* piv,
                                                      int pstride,
                                                      thrust::complex<T> LogOverlapFactor,
                                                      thrust::complex<T>* det,
                                                      int nbatch)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 64> tmp;
  int batch = blockIdx.x;
  if (batch >= nbatch)
    return;
  int t = threadIdx.x;

  tmp[t] = thrust::complex<T>(0.0);

  int const* piv_              = piv + batch * pstride;
  thrust::complex<T> const* m_ = m[batch];

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
    if (piv_[ip] == (ip + 1))
    {
      tmp[t] += log(m_[ip * lda + ip]);
    }
    else
    {
      tmp[t] += log(-m_[ip * lda + ip]);
    }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
      tmp[0] += tmp[i];
    *(det + batch) = exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<typename T>
__global__ void kernel_determinant_from_geqrf(int N,
                                              T* m,
                                              int lda,
                                              T* buff,
                                              T LogOverlapFactor,
                                              thrust::complex<T>* det)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<T, 256> tmp;
  int t = threadIdx.x;

  tmp[t] = T(0.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m[ip * lda + ip] < 0)
      buff[ip] = T(-1.0);
    else
      buff[ip] = T(1.0);
    tmp[t] += log(buff[ip] * m[ip * lda + ip]);
  }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
      tmp[0] += tmp[i];
    *det = thrust::complex<T>(exp(tmp[0] - LogOverlapFactor), 0.0);
  }
  __syncthreads();
}

template<typename T>
__global__ void kernel_determinant_from_geqrf(int N,
                                              thrust::complex<T>* m,
                                              int lda,
                                              thrust::complex<T>* buff,
                                              thrust::complex<T> LogOverlapFactor,
                                              thrust::complex<T>* det)
{
  __shared__ thrust::cuda_cub::core::uninitialized_array<thrust::complex<T>, 256> tmp;
  int t = threadIdx.x;

  tmp[t] = thrust::complex<T>(0.0);

  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m[ip * lda + ip].real() < 0)
      buff[ip] = thrust::complex<T>(-1.0);
    else
      buff[ip] = thrust::complex<T>(1.0);
    tmp[t] += thrust::log(buff[ip] * m[ip * lda + ip]);
  }
  __syncthreads();

  // not optimal but ok for now
  if (threadIdx.x == 0)
  {
    int imax = (N > blockDim.x) ? blockDim.x : N;
    for (int i = 1; i < imax; i++)
      tmp[0] += tmp[i];
    *det = thrust::exp(tmp[0] - LogOverlapFactor);
  }
  __syncthreads();
}

template<typename T>
__global__ void kernel_determinant_from_geqrf(int N, T* m, int lda, T* buff)
{
  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m[ip * lda + ip] < 0)
      buff[ip] = T(-1.0);
    else
      buff[ip] = T(1.0);
  }
}

template<typename T>
__global__ void kernel_determinant_from_geqrf(int N, thrust::complex<T>* m, int lda, thrust::complex<T>* buff)
{
  for (int ip = threadIdx.x; ip < N; ip += blockDim.x)
  {
    if (m[ip * lda + ip].real() < 0)
      buff[ip] = thrust::complex<T>(-1.0);
    else
      buff[ip] = thrust::complex<T>(1.0);
  }
}

template<typename T>
__global__ void kernel_scale_columns(int n, int m, T* A, int lda, T* scl)
{
  for (int ip = threadIdx.x; ip < n; ip += blockDim.x)
    for (int jp = threadIdx.y; jp < m; jp += blockDim.y)
      A[ip * lda + jp] *= scl[jp];
}

void determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double LogOverlapFactor, double* res)
{
  kernel_determinant_from_getrf<<<1, 256>>>(N, m, lda, piv, LogOverlapFactor, res);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void determinant_from_getrf_gpu(int N,
                                std::complex<double>* m,
                                int lda,
                                int* piv,
                                std::complex<double> LogOverlapFactor,
                                std::complex<double>* res)
{
  kernel_determinant_from_getrf<<<1, 256>>>(N, reinterpret_cast<thrust::complex<double>*>(m), lda, piv,
                                            static_cast<thrust::complex<double>>(LogOverlapFactor),
                                            reinterpret_cast<thrust::complex<double>*>(res));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void strided_determinant_from_getrf_gpu(int N,
                                        double* m,
                                        int lda,
                                        int mstride,
                                        int* piv,
                                        int pstride,
                                        double LogOverlapFactor,
                                        double* res,
                                        int nbatch)
{
  kernel_strided_determinant_from_getrf<<<nbatch, 64>>>(N, m, lda, mstride, piv, pstride, LogOverlapFactor, res,
                                                        nbatch);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void strided_determinant_from_getrf_gpu(int N,
                                        std::complex<double>* m,
                                        int lda,
                                        int mstride,
                                        int* piv,
                                        int pstride,
                                        std::complex<double> LogOverlapFactor,
                                        std::complex<double>* res,
                                        int nbatch)
{
  kernel_strided_determinant_from_getrf<<<nbatch, 64>>>(N, reinterpret_cast<thrust::complex<double>*>(m), lda, mstride,
                                                        piv, pstride,
                                                        static_cast<thrust::complex<double>>(LogOverlapFactor),
                                                        reinterpret_cast<thrust::complex<double>*>(res), nbatch);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void batched_determinant_from_getrf_gpu(int N,
                                        double** m,
                                        int lda,
                                        int* piv,
                                        int pstride,
                                        double LogOverlapFactor,
                                        double* res,
                                        int nbatch)
{
  kernel_batched_determinant_from_getrf<<<nbatch, 64>>>(N, m, lda, piv, pstride, LogOverlapFactor, res, nbatch);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void batched_determinant_from_getrf_gpu(int N,
                                        std::complex<double>** m,
                                        int lda,
                                        int* piv,
                                        int pstride,
                                        std::complex<double> LogOverlapFactor,
                                        std::complex<double>* res,
                                        int nbatch)
{
  kernel_batched_determinant_from_getrf<<<nbatch, 64>>>(N, reinterpret_cast<thrust::complex<double>**>(m), lda, piv,
                                                        pstride, static_cast<thrust::complex<double>>(LogOverlapFactor),
                                                        reinterpret_cast<thrust::complex<double>*>(res), nbatch);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

std::complex<double> determinant_from_geqrf_gpu(int N, double* m, int lda, double* buff, double LogOverlapFactor)
{
  thrust::device_ptr<thrust::complex<double>> d_ptr = thrust::device_malloc<thrust::complex<double>>(1);
  kernel_determinant_from_geqrf<<<1, 256>>>(N, m, lda, buff, LogOverlapFactor, thrust::raw_pointer_cast(d_ptr));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  std::complex<double> res;
  qmc_cuda::cuda_check(cudaMemcpy(std::addressof(res), thrust::raw_pointer_cast(d_ptr), sizeof(std::complex<double>),
                                  cudaMemcpyDeviceToHost));
  thrust::device_free(d_ptr);
  return res;
}

std::complex<double> determinant_from_geqrf_gpu(int N,
                                                std::complex<double>* m,
                                                int lda,
                                                std::complex<double>* buff,
                                                std::complex<double> LogOverlapFactor)
{
  thrust::device_ptr<thrust::complex<double>> d_ptr = thrust::device_malloc<thrust::complex<double>>(1);
  kernel_determinant_from_geqrf<<<1, 256>>>(N, reinterpret_cast<thrust::complex<double>*>(m), lda,
                                            reinterpret_cast<thrust::complex<double>*>(buff),
                                            static_cast<thrust::complex<double>>(LogOverlapFactor),
                                            thrust::raw_pointer_cast(d_ptr));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  std::complex<double> res;
  qmc_cuda::cuda_check(cudaMemcpy(std::addressof(res), thrust::raw_pointer_cast(d_ptr), sizeof(std::complex<double>),
                                  cudaMemcpyDeviceToHost));
  thrust::device_free(d_ptr);
  return res;
}

void determinant_from_geqrf_gpu(int N, double* m, int lda, double* buff)
{
  kernel_determinant_from_geqrf<<<1, 256>>>(N, m, lda, buff);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

void determinant_from_geqrf_gpu(int N, std::complex<double>* m, int lda, std::complex<double>* buff)
{
  kernel_determinant_from_geqrf<<<1, 256>>>(N, reinterpret_cast<thrust::complex<double>*>(m), lda,
                                            reinterpret_cast<thrust::complex<double>*>(buff));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}


double determinant_from_getrf_gpu(int N, double* m, int lda, int* piv, double LogOverlapFactor)
{
  thrust::device_ptr<double> d_ptr = thrust::device_malloc<double>(1);
  kernel_determinant_from_getrf<<<1, 256>>>(N, m, lda, piv, LogOverlapFactor, thrust::raw_pointer_cast(d_ptr));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  double res = *d_ptr;
  thrust::device_free(d_ptr);
  return res;
}

std::complex<double> determinant_from_getrf_gpu(int N,
                                                std::complex<double>* m,
                                                int lda,
                                                int* piv,
                                                std::complex<double> LogOverlapFactor)
{
  thrust::device_ptr<thrust::complex<double>> d_ptr = thrust::device_malloc<thrust::complex<double>>(1);
  kernel_determinant_from_getrf<<<1, 256>>>(N, reinterpret_cast<thrust::complex<double>*>(m), lda, piv,
                                            static_cast<thrust::complex<double>>(LogOverlapFactor),
                                            thrust::raw_pointer_cast(d_ptr));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
  std::complex<double> res;
  qmc_cuda::cuda_check(cudaMemcpy(std::addressof(res), thrust::raw_pointer_cast(d_ptr), sizeof(std::complex<double>),
                                  cudaMemcpyDeviceToHost));
  thrust::device_free(d_ptr);
  return res;
}

void scale_columns(int n, int m, double* A, int lda, double* scl)
{
  int xblock_dim = 32;
  dim3 block_dim(xblock_dim, xblock_dim, 1);
  kernel_scale_columns<<<1, block_dim>>>(n, m, A, lda, scl);
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}
void scale_columns(int n, int m, std::complex<double>* A, int lda, std::complex<double>* scl)
{
  int xblock_dim = 32;
  dim3 block_dim(xblock_dim, xblock_dim, 1);
  kernel_scale_columns<<<1, block_dim>>>(n, m, reinterpret_cast<thrust::complex<double>*>(A), lda,
                                         reinterpret_cast<thrust::complex<double>*>(scl));
  qmc_cuda::cuda_check(cudaGetLastError());
  qmc_cuda::cuda_check(cudaDeviceSynchronize());
}

} // namespace kernels
