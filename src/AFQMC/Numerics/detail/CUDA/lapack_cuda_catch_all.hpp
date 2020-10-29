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

#ifndef AFQMC_LAPACK_CUDA_CATCH_ALL_HPP
#define AFQMC_LAPACK_CUDA_CATCH_ALL_HPP

#include <cassert>
#include "AFQMC/Utilities/type_conversion.hpp"
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/detail/CPU/lapack_cpu.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/cusolver_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/Kernels/setIdentity.cuh"

namespace device
{
using qmcplusplus::afqmc::remove_complex;

// hevr
template<typename T, class ptr, class ptrR, class ptrI>
inline static void hevr(char JOBZ,
                        char RANGE,
                        char UPLO,
                        int N,
                        ptr A,
                        int LDA,
                        T VL,
                        T VU,
                        int IL,
                        int IU,
                        T ABSTOL,
                        int& M,
                        ptrR W,
                        ptr Z,
                        int LDZ,
                        ptrI ISUPPZ,
                        ptr WORK,
                        int& LWORK,
                        ptrR RWORK,
                        int& LRWORK,
                        ptrI IWORK,
                        int& LIWORK,
                        int& INFO)
{
  throw std::runtime_error("Error: Calling qmc_cuda::hevr catch all.");
}

// getrf
template<class ptr, class ptrW, class ptrI>
inline static void getrf(const int n, const int m, ptr const a, int lda, ptrI piv, int& st, ptrW work)
{
  throw std::runtime_error("Error: Calling qmc_cuda::getrf catch all.");
}

// getrfBatched
template<class ptr, class ptrI1, class ptrI2>
inline static void getrfBatched(const int n, ptr* a, int lda, ptrI1 piv, ptrI2 info, int batchSize)
{
  throw std::runtime_error("Error: Calling qmc_cuda::getrfBatched catch all.");
}

// getri: will fail if not called correctly, but removing checks on ptrI and ptrW for now
template<class ptr, class ptrI, class ptrW>
inline static void getri(int n, ptr a, int n0, ptrI piv, ptrW work, int n1, int& status)
{
  throw std::runtime_error("Error: Calling qmc_cuda::getri catch all.");
}

// getriBatched
template<class ptr, class ptrR, class ptrI1, class ptrI2>
inline static void getriBatched(int n, ptr* a, int lda, ptrI1 piv, ptrR* c, int lwork, ptrI2 info, int batchSize)
{
  throw std::runtime_error("Error: Calling qmc_cuda::getriBatched catch all.");
}

template<class ptrA, class ptrB, class ptrC>
inline static void geqrf(int M, int N, ptrA A, const int LDA, ptrB TAU, ptrC WORK, int LWORK, int& INFO)
{
  throw std::runtime_error("Error: Calling qmc_cuda::geqrf catch all.");
}

template<class ptrA, class ptrB, class ptrC>
void static gqr(int M, int N, int K, ptrA A, const int LDA, ptrB TAU, ptrC WORK, int LWORK, int& INFO)
{
  throw std::runtime_error("Error: Calling qmc_cuda::gqr catch all.");
}

template<class ptrA, class ptrB, class ptrC>
inline static void gelqf(int M, int N, ptrA A, const int LDA, ptrB TAU, ptrC WORK, int LWORK, int& INFO)
{
  throw std::runtime_error("Error: Calling qmc_cuda::gelqf catch all.");
}

template<class ptrA, class ptrB, class ptrC>
void static glq(int M, int N, int K, ptrA A, const int LDA, ptrB TAU, ptrC WORK, int LWORK, int& INFO)
{
  throw std::runtime_error("Error: Calling qmc_cuda::glq catch all.");
}

} // namespace device

#endif
