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

#ifndef AFQMC_NUMERICS_HELPERS_HPP
#define AFQMC_NUMERICS_HELPERS_HPP

#include <cassert>
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/device_kernels.hpp"
#endif

namespace ma
{
template<class T>
inline T determinant_from_getrf(int n, T* M, int lda, int* pivot, T LogOverlapFactor)
{
  T res(0.0);
  T sg(1.0);
  for (int i = 0, ip = 1; i != n; i++, ip++)
  {
    if (real(M[i * lda + i]) < 0.0)
    {
      res += std::log(-static_cast<T>(M[i * lda + i]));
      sg *= -1.0;
    }
    else
      res += std::log(static_cast<T>(M[i * lda + i]));
    if (pivot[i] != ip)
      sg *= -1.0;
  }
  return sg * std::exp(res - LogOverlapFactor);
}

template<class T>
inline void determinant_from_getrf(int n, T* M, int lda, int* pivot, T LogOverlapFactor, T* res)
{
  *res = T(0.0);
  T sg(1.0);
  for (int i = 0, ip = 1; i != n; i++, ip++)
  {
    if (real(M[i * lda + i]) < 0.0)
    {
      *res += std::log(-static_cast<T>(M[i * lda + i]));
      sg *= -1.0;
    }
    else
      *res += std::log(static_cast<T>(M[i * lda + i]));
    if (pivot[i] != ip)
      sg *= -1.0;
  }
  *res = sg * std::exp(*res - LogOverlapFactor);
}

template<class T>
inline void strided_determinant_from_getrf(int n,
                                           T* M,
                                           int lda,
                                           int Mstride,
                                           int* pivot,
                                           int pstride,
                                           T LogOverlapFactor,
                                           T* res,
                                           int nbatch)
{
  for (int b = 0; b < nbatch; b++)
    determinant_from_getrf(n, M + b * Mstride, lda, pivot + b * pstride, LogOverlapFactor, res + b);
}

template<class T>
inline void batched_determinant_from_getrf(int n,
                                           T** M,
                                           int lda,
                                           int* pivot,
                                           int pstride,
                                           T LogOverlapFactor,
                                           T* res,
                                           int nbatch)
{
  for (int b = 0; b < nbatch; b++)
    determinant_from_getrf(n, M[b], lda, pivot + b * pstride, LogOverlapFactor, res + b);
}

template<class T>
T determinant_from_geqrf(int n, T* M, int lda, T* buff, T LogOverlapFactor)
{
  T res(0.0);
  for (int i = 0; i < n; i++)
  {
    if (real(M[i * lda + i]) < 0.0)
      buff[i] = T(-1.0);
    else
      buff[i] = T(1.0);
    res += std::log(buff[i] * M[i * lda + i]);
  }
  return std::exp(res - LogOverlapFactor);
}

// specializations for complex
template<class T>
inline std::complex<T> determinant_from_getrf(int n,
                                              std::complex<T>* M,
                                              int lda,
                                              int* pivot,
                                              std::complex<T> LogOverlapFactor)
{
  std::complex<T> res(0.0, 0.0);
  for (int i = 0, ip = 1; i != n; i++, ip++)
  {
    if (pivot[i] == ip)
    {
      res += std::log(+static_cast<std::complex<T>>(M[i * lda + i]));
    }
    else
    {
      res += std::log(-static_cast<std::complex<T>>(M[i * lda + i]));
    }
  }
  return std::exp(res - LogOverlapFactor);
}

template<class T>
inline void determinant_from_getrf(int n,
                                   std::complex<T>* M,
                                   int lda,
                                   int* pivot,
                                   std::complex<T> LogOverlapFactor,
                                   std::complex<T>* res)
{
  *res = std::complex<T>(0.0, 0.0);
  for (int i = 0, ip = 1; i != n; i++, ip++)
  {
    if (pivot[i] == ip)
    {
      *res += std::log(+static_cast<std::complex<T>>(M[i * lda + i]));
    }
    else
    {
      *res += std::log(-static_cast<std::complex<T>>(M[i * lda + i]));
    }
  }
  *res = std::exp(*res - LogOverlapFactor);
}

template<class T>
inline void determinant_from_geqrf(int n, T* M, int lda, T* buff)
{
  for (int i = 0; i < n; i++)
  {
    if (real(M[i * lda + i]) < 0)
      buff[i] = T(-1.0);
    else
      buff[i] = T(1.0);
  }
}

template<class T>
inline void scale_columns(int n, int m, T* A, int lda, T* scl)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      A[i * lda + j] *= scl[j];
}

} // namespace ma

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
namespace device
{
// using thrust for now to avoid kernels!!!
template<class T>
inline void determinant_from_getrf(int n,
                                   device_pointer<T> A,
                                   int lda,
                                   device_pointer<int> piv,
                                   T LogOverlapFactor,
                                   T* res)
{
  kernels::determinant_from_getrf_gpu(n, to_address(A), lda, to_address(piv), LogOverlapFactor, res);
}

template<class T>
inline void strided_determinant_from_getrf(int n,
                                           device_pointer<T> A,
                                           int lda,
                                           int Mstride,
                                           device_pointer<int> piv,
                                           int pstride,
                                           T LogOverlapFactor,
                                           T* res,
                                           int nbatch)
{
  kernels::strided_determinant_from_getrf_gpu(n, to_address(A), lda, Mstride, to_address(piv), pstride,
                                              LogOverlapFactor, res, nbatch);
}

template<class T>
inline void batched_determinant_from_getrf(int n,
                                           device_pointer<T>* A,
                                           int lda,
                                           device_pointer<int> piv,
                                           int pstride,
                                           T LogOverlapFactor,
                                           T* res,
                                           int nbatch)
{
  T** A_h = new T*[nbatch];
  for (int i = 0; i < nbatch; i++)
    A_h[i] = to_address(A[i]);
  T** A_d;
  arch::malloc((void**)&A_d, nbatch * sizeof(*A_h));
  arch::memcopy(A_d, A_h, nbatch * sizeof(*A_h), arch::memcopyH2D);
  kernels::batched_determinant_from_getrf_gpu(n, A_d, lda, to_address(piv), pstride, LogOverlapFactor, res, nbatch);
  arch::free(A_d);
  delete[] A_h;
}

template<class T>
inline T determinant_from_getrf(int n, device_pointer<T> A, int lda, device_pointer<int> piv, T LogOverlapFactor)
{
  return kernels::determinant_from_getrf_gpu(n, to_address(A), lda, to_address(piv), LogOverlapFactor);
}

template<class T>
T determinant_from_geqrf(int n, device_pointer<T> M, int lda, device_pointer<T> buff, T LogOverlapFactor)
{
  return kernels::determinant_from_geqrf_gpu(n, to_address(M), lda, to_address(buff), LogOverlapFactor);
}

template<class T>
inline void determinant_from_geqrf(int n, device_pointer<T> M, int lda, device_pointer<T> buff)
{
  kernels::determinant_from_geqrf_gpu(n, to_address(M), lda, to_address(buff));
}

template<class T>
inline void scale_columns(int n, int m, device_pointer<T> A, int lda, device_pointer<T> scl)
{
  kernels::scale_columns(n, m, to_address(A), lda, to_address(scl));
}

template<class ptrA, class ptrB>
inline void scale_columns(int n, int m, ptrA A, int lda, ptrB scl)
{
  print_stacktrace;
  throw std::runtime_error("Error: Calling device::scale_columns atch all.");
}

} // namespace device
#endif

#endif
