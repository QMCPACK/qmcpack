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
#include <hip/hip_runtime.h>
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include "AFQMC/Numerics/detail/HIP/Kernels/hip_settings.h"
#include "AFQMC/Numerics/detail/HIP/hip_kernel_utils.h"

namespace kernels
{
// very sloppy, needs improvement!!!!
// v1[nkpts][nchol_max][nwalk]
// vb[2*nchol_tot][nwalk]
// ncholpQ0 includes factor of 2
template<typename T, typename T2>
__global__ void kernel_vbias_from_v1(int nwalk,
                                     int nkpts,
                                     int nchol_max,
                                     int* Qsym,
                                     int* kminus,
                                     int* ncholpQ,
                                     int* ncholpQ0,
                                     thrust::complex<T2> const alpha,
                                     thrust::complex<T> const* v1,
                                     thrust::complex<T2>* vb)
{
  // This is inefficient now, eventually map to assigned Qs and only launch based on this
  // essentially only need a map to the assigned Q's and use blockIdx.x to index into this array
  int Q = blockIdx.x;
  if (Q >= nkpts || blockIdx.y > 1)
    return;
  if (Qsym[Q] < 0)
    return;
  int Qm  = kminus[Q];
  int nc0 = ncholpQ0[Q];
  int nc  = ncholpQ[Q];
  int ncm = ncholpQ[Qm];
  // now redefine Qm based on Qsym
  if (Qsym[Q] > 0)
    Qm = nkpts + Qsym[Q] - 1;

  if (blockIdx.y == 0)
  {
    // v+
    thrust::complex<T2>* vb_(vb + nc0 * nwalk);
    thrust::complex<T> const* v1_(v1 + Q * nchol_max * nwalk);
    thrust::complex<T> const* v2_(v1 + Qm * nchol_max * nwalk);
    // v+ = a*(v[Q]+v[-Q])
    if (threadIdx.x < nc && threadIdx.y < nwalk)
    {
      for (int n = threadIdx.x; n < nc; n += blockDim.x)
        for (int w = threadIdx.y; w < nwalk; w += blockDim.y)
          vb_[n * nwalk + w] += alpha * static_cast<thrust::complex<T2>>(v1_[n * nwalk + w]);
    }
    if (threadIdx.x < ncm && threadIdx.y < nwalk)
    {
      for (int n = threadIdx.x; n < ncm; n += blockDim.x)
        for (int w = threadIdx.y; w < nwalk; w += blockDim.y)
          vb_[n * nwalk + w] += alpha * static_cast<thrust::complex<T2>>(v2_[n * nwalk + w]);
    }
  }
  else if (blockIdx.y == 1)
  {
    // v-
    thrust::complex<T2>* vb_(vb + (nc0 + nc) * nwalk);
    thrust::complex<T> const* v1_(v1 + Q * nchol_max * nwalk);
    thrust::complex<T> const* v2_(v1 + Qm * nchol_max * nwalk);
    // v- = -a*i*(v[Q]-v[-Q])
    thrust::complex<T2> ialpha(alpha * thrust::complex<T2>(0.0, 1.0));
    if (threadIdx.x < nc && threadIdx.y < nwalk)
    {
      for (int n = threadIdx.x; n < nc; n += blockDim.x)
        for (int w = threadIdx.y; w < nwalk; w += blockDim.y)
          vb_[n * nwalk + w] -= ialpha * static_cast<thrust::complex<T2>>(v1_[n * nwalk + w]);
    }
    if (threadIdx.x < ncm && threadIdx.y < nwalk)
    {
      for (int n = threadIdx.x; n < ncm; n += blockDim.x)
        for (int w = threadIdx.y; w < nwalk; w += blockDim.y)
          vb_[n * nwalk + w] += ialpha * static_cast<thrust::complex<T2>>(v2_[n * nwalk + w]);
    }
  }
}

void vbias_from_v1(int nwalk,
                   int nkpts,
                   int nchol_max,
                   int* Qsym,
                   int* kminus,
                   int* ncholpQ,
                   int* ncholpQ0,
                   std::complex<double> const alpha,
                   std::complex<double> const* v1,
                   std::complex<double>* vb)
{
  int xblock_dim = 32;
  int yblock_dim = std::min(nwalk, 16);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, 2, 1);
  hipLaunchKernelGGL(kernel_vbias_from_v1, dim3(grid_dim), dim3(block_dim), 0, 0, nwalk, nkpts, nchol_max, Qsym, kminus,
                     ncholpQ, ncholpQ0, static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<double> const*>(v1),
                     reinterpret_cast<thrust::complex<double>*>(vb));
  qmc_hip::hip_kernel_check(hipGetLastError(), "vbias_from_v1");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "vbias_from_v1");
}

void vbias_from_v1(int nwalk,
                   int nkpts,
                   int nchol_max,
                   int* Qsym,
                   int* kminus,
                   int* ncholpQ,
                   int* ncholpQ0,
                   std::complex<float> const alpha,
                   std::complex<float> const* v1,
                   std::complex<float>* vb)
{
  int xblock_dim = 32;
  int yblock_dim = std::min(nwalk, 16);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, 2, 1);
  hipLaunchKernelGGL(kernel_vbias_from_v1, dim3(grid_dim), dim3(block_dim), 0, 0, nwalk, nkpts, nchol_max, Qsym, kminus,
                     ncholpQ, ncholpQ0, static_cast<thrust::complex<float> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(v1),
                     reinterpret_cast<thrust::complex<float>*>(vb));
  qmc_hip::hip_kernel_check(hipGetLastError(), "vbias_from_v1");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "vbias_from_v1");
}

void vbias_from_v1(int nwalk,
                   int nkpts,
                   int nchol_max,
                   int* Qsym,
                   int* kminus,
                   int* ncholpQ,
                   int* ncholpQ0,
                   std::complex<double> const alpha,
                   std::complex<float> const* v1,
                   std::complex<double>* vb)
{
  int xblock_dim = 32;
  int yblock_dim = std::min(nwalk, 16);
  dim3 block_dim(xblock_dim, yblock_dim, 1);
  dim3 grid_dim(nkpts, 2, 1);
  hipLaunchKernelGGL(kernel_vbias_from_v1, dim3(grid_dim), dim3(block_dim), 0, 0, nwalk, nkpts, nchol_max, Qsym, kminus,
                     ncholpQ, ncholpQ0, static_cast<thrust::complex<double> const>(alpha),
                     reinterpret_cast<thrust::complex<float> const*>(v1),
                     reinterpret_cast<thrust::complex<double>*>(vb));
  qmc_hip::hip_kernel_check(hipGetLastError(), "vbias_from_v1");
  qmc_hip::hip_kernel_check(hipDeviceSynchronize(), "vbias_from_v1");
}

} // namespace kernels
