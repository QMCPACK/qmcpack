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
#include <CL/sycl.hpp>

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
} // namespace syclBLAS

} // namespace qmcplusplus
#endif // QMCPLUSPLUS_OMPBLAS_H
