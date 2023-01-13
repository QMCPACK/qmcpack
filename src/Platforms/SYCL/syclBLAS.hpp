//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SYCL_BLAS_H
#define QMCPLUSPLUS_SYCL_BLAS_H

#include <complex>
#include <sycl/sycl.hpp>

namespace qmcplusplus
{
namespace syclBLAS
{
using syclBLAS_int    = std::int64_t;
using syclBLAS_status = sycl::event;
using syclBLAS_handle = sycl::queue;

template<typename T>
sycl::event gemv(sycl::queue& handle,
                 const char trans,
                 const int m,
                 const int n,
                 const T alpha,
                 const T* const A,
                 const int lda,
                 const T* const x,
                 const int incx,
                 const T beta,
                 T* const y,
                 const int incy,
                 const std::vector<sycl::event>& events = {});

template<typename T>
sycl::event gemm(sycl::queue& handle,
                 const char tA,
                 const char tB,
                 const int m,
                 const int n,
                 const int k,
                 const T alpha,
                 const T* const A,
                 const int lda,
                 const T* const B,
                 const int ldb,
                 const T beta,
                 T* const C,
                 const int ldc,
                 const std::vector<sycl::event>& events = {});

template<typename T1, typename T2>
sycl::event transpose(sycl::queue& q,
                      const T1* in,
                      int m,
                      int lda,
                      T2* out,
                      int n,
                      int ldb,
                      const std::vector<sycl::event>& events = {});

template<typename T1, typename T2>
sycl::event copy_n(sycl::queue& aq,
                   const T1* VA,
                   size_t array_size,
                   T2* VC,
                   const std::vector<sycl::event>& events = {});

} // namespace syclBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_OMPBLAS_H
