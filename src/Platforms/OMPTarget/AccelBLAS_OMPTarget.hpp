//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
// File refactored from: MatrixDelayedUpdateCUDA.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ACCELBLAS_OMPTARGET_H
#define QMCPLUSPLUS_ACCELBLAS_OMPTARGET_H

#include "AccelBLAS.hpp"
#include "QueueOMPTarget.hpp"
#include "ompBLAS.hpp"

namespace qmcplusplus
{
namespace compute
{
template<>
class BLASHandle<PlatformKind::OMPTARGET>
{
public:
  ompBLAS::ompBLAS_handle h_ompblas;

  BLASHandle(Queue<PlatformKind::OMPTARGET>& queue) : h_ompblas(0) { }
};

namespace BLAS
{
template<>
inline void gemm<PlatformKind::OMPTARGET, float>(BLASHandle<PlatformKind::OMPTARGET>& handle,
                                            const char transa,
                                            const char transb,
                                            int m,
                                            int n,
                                            int k,
                                            const float alpha,
                                            const float* A,
                                            int lda,
                                            const float* B,
                                            int ldb,
                                            const float beta,
                                            float* C,
                                            int ldc)
{
  if(ompBLAS::gemm(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) != 0)
    throw std::runtime_error("ompBLAS::gemm<float> failed!");
}

template<>
inline void gemm<PlatformKind::OMPTARGET, double>(BLASHandle<PlatformKind::OMPTARGET>& handle,
                                             const char transa,
                                             const char transb,
                                             int m,
                                             int n,
                                             int k,
                                             const double alpha,
                                             const double* A,
                                             int lda,
                                             const double* B,
                                             int ldb,
                                             const double beta,
                                             double* C,
                                             int ldc)
{
  if(ompBLAS::gemm(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) != 0)
    throw std::runtime_error("ompBLAS::gemm<double> failed!");
}

template<>
inline void gemm<PlatformKind::OMPTARGET, std::complex<float>>(BLASHandle<PlatformKind::OMPTARGET>& handle,
                                                          const char transa,
                                                          const char transb,
                                                          int m,
                                                          int n,
                                                          int k,
                                                          const std::complex<float> alpha,
                                                          const std::complex<float>* A,
                                                          int lda,
                                                          const std::complex<float>* B,
                                                          int ldb,
                                                          const std::complex<float> beta,
                                                          std::complex<float>* C,
                                                          int ldc)
{
  if(ompBLAS::gemm(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) != 0)
    throw std::runtime_error("ompBLAS::gemm<std::complex<float>> failed!");
}

template<>
inline void gemm<PlatformKind::OMPTARGET, std::complex<double>>(BLASHandle<PlatformKind::OMPTARGET>& handle,
                                                           const char transa,
                                                           const char transb,
                                                           int m,
                                                           int n,
                                                           int k,
                                                           const std::complex<double> alpha,
                                                           const std::complex<double>* A,
                                                           int lda,
                                                           const std::complex<double>* B,
                                                           int ldb,
                                                           const std::complex<double> beta,
                                                           std::complex<double>* C,
                                                           int ldc)
{
  if(ompBLAS::gemm(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) != 0)
    throw std::runtime_error("ompBLAS::gemm<std::complex<double>> failed!");
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
