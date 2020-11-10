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
// Tab [nbatch][nwalk][nocc][nocc][nchol]
// Klr [nwalk][2*nchol_tot]
template<typename T>
__global__ void kernel_batched_Tab_to_Klr(int nterms,
                                          int nwalk,
                                          int nocc,
                                          int nchol_max,
                                          int nchol_tot,
                                          int ncholQ,
                                          int ncholQ0,
                                          int* kdiag,
                                          thrust::complex<T> const* Tab,
                                          thrust::complex<T>* Kl,
                                          thrust::complex<T>* Kr)
{
  int w = blockIdx.x;
  if (blockIdx.y == 0)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      if (w < nwalk)
      {
        for (int a = 0; a < nocc; a++)
        {
          thrust::complex<T> const* Tba_(Tab + batch * nwalk * nocc * nocc * nchol_max +
                                         ((w * nocc + a) * nocc + a) * nchol_max);
          thrust::complex<T>* Kr_(Kr + w * nchol_tot + ncholQ0);
          int c = threadIdx.x;
          while (c < ncholQ)
          {
            Kr_[c] += Tba_[c];
            c += blockDim.x;
          }
        }
      }
    }
  }
  else if (blockIdx.y == 1)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      if (w < nwalk)
      {
        for (int a = 0; a < nocc; a++)
        {
          thrust::complex<T> const* Tab_(Tab + (batch + 1) * nwalk * nocc * nocc * nchol_max +
                                         ((w * nocc + a) * nocc + a) * nchol_max);
          thrust::complex<T>* Kl_(Kl + w * nchol_tot + ncholQ0);
          int c = threadIdx.x;
          while (c < ncholQ)
          {
            Kl_[c] += Tab_[c];
            c += blockDim.x;
          }
        }
      }
    }
  }
}

template<typename T>
__global__ void kernel_batched_Tanb_to_Klr(int nterms,
                                           int nwalk,
                                           int nocc,
                                           int nchol_max,
                                           int nchol_tot,
                                           int ncholQ,
                                           int ncholQ0,
                                           int* kdiag,
                                           thrust::complex<T> const* Tab,
                                           thrust::complex<T>* Kl,
                                           thrust::complex<T>* Kr)
{
  int w = blockIdx.x;
  if (blockIdx.y == 0)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      if (w < nwalk)
      {
        for (int a = 0; a < nocc; a++)
        {
          thrust::complex<T> const* Tba_(Tab + batch * nwalk * nocc * nocc * nchol_max +
                                         ((w * nocc + a) * nocc) * nchol_max + a);
          thrust::complex<T>* Kr_(Kr + w * nchol_tot + ncholQ0);
          int c = threadIdx.x;
          while (c < ncholQ)
          {
            Kr_[c] += Tba_[c * nocc];
            c += blockDim.x;
          }
        }
      }
    }
  }
  else if (blockIdx.y == 1)
  {
    for (int k = 0; k < nterms; k++)
    {
      int batch = kdiag[k];
      if (w < nwalk)
      {
        for (int a = 0; a < nocc; a++)
        {
          thrust::complex<T> const* Tab_(Tab + (batch + 1) * nwalk * nocc * nocc * nchol_max +
                                         ((w * nocc + a) * nocc) * nchol_max + a);
          thrust::complex<T>* Kl_(Kl + w * nchol_tot + ncholQ0);
          int c = threadIdx.x;
          while (c < ncholQ)
          {
            Kl_[c] += Tab_[c * nocc];
            c += blockDim.x;
          }
        }
      }
    }
  }
}

void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        int* kdiag,
                        std::complex<double> const* Tab,
                        std::complex<double>* Kl,
                        std::complex<double>* Kr)
{
  dim3 grid_dim(nwalk, 2, 1);
  int nthr = std::min(256, ncholQ); // is this needed?
  kernel_batched_Tab_to_Klr<<<grid_dim, nthr>>>(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, kdiag,
                                                reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                reinterpret_cast<thrust::complex<double>*>(Kl),
                                                reinterpret_cast<thrust::complex<double>*>(Kr));
  qmc_cuda::cuda_check(cudaGetLastError(), "batched_Tab_to_Klr");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "batched_Tab_to_Klr");
}

void batched_Tab_to_Klr(int nterms,
                        int nwalk,
                        int nocc,
                        int nchol_max,
                        int nchol_tot,
                        int ncholQ,
                        int ncholQ0,
                        int* kdiag,
                        std::complex<float> const* Tab,
                        std::complex<float>* Kl,
                        std::complex<float>* Kr)
{
  dim3 grid_dim(nwalk, 2, 1);
  int nthr = std::min(256, ncholQ); // is this needed?
  kernel_batched_Tab_to_Klr<<<grid_dim, nthr>>>(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, kdiag,
                                                reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                reinterpret_cast<thrust::complex<float>*>(Kl),
                                                reinterpret_cast<thrust::complex<float>*>(Kr));
  qmc_cuda::cuda_check(cudaGetLastError(), "batched_Tab_to_Klr");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "batched_Tab_to_Klr");
}

void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         int* kdiag,
                         std::complex<double> const* Tab,
                         std::complex<double>* Kl,
                         std::complex<double>* Kr)
{
  dim3 grid_dim(nwalk, 2, 1);
  int nthr = std::min(256, ncholQ); // is this needed?
  kernel_batched_Tanb_to_Klr<<<grid_dim, nthr>>>(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, kdiag,
                                                 reinterpret_cast<thrust::complex<double> const*>(Tab),
                                                 reinterpret_cast<thrust::complex<double>*>(Kl),
                                                 reinterpret_cast<thrust::complex<double>*>(Kr));
  qmc_cuda::cuda_check(cudaGetLastError(), "batched_Tanb_to_Klr");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "batched_Tanb_to_Klr");
}

void batched_Tanb_to_Klr(int nterms,
                         int nwalk,
                         int nocc,
                         int nchol_max,
                         int nchol_tot,
                         int ncholQ,
                         int ncholQ0,
                         int* kdiag,
                         std::complex<float> const* Tab,
                         std::complex<float>* Kl,
                         std::complex<float>* Kr)
{
  dim3 grid_dim(nwalk, 2, 1);
  int nthr = std::min(256, ncholQ); // is this needed?
  kernel_batched_Tanb_to_Klr<<<grid_dim, nthr>>>(nterms, nwalk, nocc, nchol_max, nchol_tot, ncholQ, ncholQ0, kdiag,
                                                 reinterpret_cast<thrust::complex<float> const*>(Tab),
                                                 reinterpret_cast<thrust::complex<float>*>(Kl),
                                                 reinterpret_cast<thrust::complex<float>*>(Kr));
  qmc_cuda::cuda_check(cudaGetLastError(), "batched_Tanb_to_Klr");
  qmc_cuda::cuda_check(cudaDeviceSynchronize(), "batched_Tanb_to_Klr");
}

} // namespace kernels
