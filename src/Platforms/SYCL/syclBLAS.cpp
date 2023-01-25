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


#include "syclBLAS.hpp"
#include "oneapi/mkl/blas.hpp"

namespace qmcplusplus
{
namespace syclBLAS
{
inline oneapi::mkl::transpose convertTransEnum(char trans)
{
  return trans == 'T' ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::nontrans;
}

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
                 const std::vector<sycl::event>& events)
{
  return oneapi::mkl::blas::gemv(handle, convertTransEnum(trans), m, n, alpha, A, lda, x, incx, beta, y, incy, events);
}

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const double alpha,
                          const double* const A,
                          const int lda,
                          const double* const x,
                          const int incx,
                          const double beta,
                          double* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const float alpha,
                          const float* const A,
                          const int lda,
                          const float* const x,
                          const int incx,
                          const float beta,
                          float* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const std::complex<double> alpha,
                          const std::complex<double>* const A,
                          const int lda,
                          const std::complex<double>* const x,
                          const int incx,
                          const std::complex<double> beta,
                          std::complex<double>* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template sycl::event gemv(sycl::queue& handle,
                          const char trans,
                          const int m,
                          const int n,
                          const std::complex<float> alpha,
                          const std::complex<float>* const A,
                          const int lda,
                          const std::complex<float>* const x,
                          const int incx,
                          const std::complex<float> beta,
                          std::complex<float>* const y,
                          const int incy,
                          const std::vector<sycl::event>& events);

template<typename T>
sycl::event gemm(sycl::queue& handle,
                 const char tA,
                 const char tB,
                 const int m,
                 const int n,
                 const int k,
                 const T alpha,
                 const T* A,
                 const int lda,
                 const T* B,
                 const int ldb,
                 const T beta,
                 T* C,
                 const int ldc,
                 const std::vector<sycl::event>& events)
{
  return oneapi::mkl::blas::gemm(handle, convertTransEnum(tA), convertTransEnum(tB), m, n, k, alpha, A, lda, B, ldb,
                                 beta, C, ldc, events);
}


template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const float alpha,
                          const float* const A,
                          const int lda,
                          const float* const B,
                          const int ldb,
                          const float beta,
                          float* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const double alpha,
                          const double* const A,
                          const int lda,
                          const double* const B,
                          const int ldb,
                          const double beta,
                          double* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<float> alpha,
                          const std::complex<float>* const A,
                          const int lda,
                          const std::complex<float>* const B,
                          const int ldb,
                          const std::complex<float> beta,
                          std::complex<float>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

template sycl::event gemm(sycl::queue& handle,
                          const char tA,
                          const char tB,
                          const int m,
                          const int n,
                          const int k,
                          const std::complex<double> alpha,
                          const std::complex<double>* const A,
                          const int lda,
                          const std::complex<double>* const B,
                          const int ldb,
                          const std::complex<double> beta,
                          std::complex<double>* const C,
                          const int ldc,
                          const std::vector<sycl::event>& events);

//transpose
template<typename T1, typename T2>
sycl::event transpose(sycl::queue& q,
                      const T1* restrict in,
                      int m,
                      int lda,
                      T2* restrict out,
                      int n,
                      int ldb,
                      const std::vector<sycl::event>& events)
{
  constexpr size_t tile_size = 16;
  const size_t m_max         = ((m + tile_size - 1) / tile_size) * tile_size;
  const size_t n_max         = ((n + tile_size - 1) / tile_size) * tile_size;

  return q.submit([&](sycl::handler& cgh) {
    cgh.depends_on(events);
    sycl::local_accessor<T2, 2> tile(sycl::range<2>(tile_size, tile_size + 1), cgh);

    cgh.parallel_for(sycl::nd_range<2>{{m_max, n_max}, {tile_size, tile_size}}, [=](sycl::nd_item<2> item) {
      unsigned x   = item.get_global_id(1);
      unsigned y   = item.get_global_id(0);
      unsigned xth = item.get_local_id(1);
      unsigned yth = item.get_local_id(0);

      if (x < n && y < m)
        tile[yth][xth] = in[(y)*lda + x];
      item.barrier(sycl::access::fence_space::local_space);

      x = item.get_group(0) * tile_size + xth;
      y = item.get_group(1) * tile_size + yth;
      if (x < m && y < n)
        out[(y)*ldb + x] = tile[xth][yth];
    });
  });
}

template sycl::event transpose(sycl::queue& q,
                               const float* restrict in,
                               int m,
                               int lda,
                               double* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const double* restrict in,
                               int m,
                               int lda,
                               double* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const std::complex<float>* restrict in,
                               int m,
                               int lda,
                               std::complex<double>* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

template sycl::event transpose(sycl::queue& q,
                               const std::complex<double>* restrict in,
                               int m,
                               int lda,
                               std::complex<double>* restrict out,
                               int n,
                               int ldb,
                               const std::vector<sycl::event>& events);

//copy_n for mixed precision
template<typename T1, typename T2>
sycl::event copy_n(sycl::queue& aq,
                   const T1* restrict VA,
                   size_t array_size,
                   T2* restrict VC,
                   const std::vector<sycl::event>& events)
{
  constexpr size_t tile_size = 64;
  const size_t a_max         = ((array_size + tile_size - 1) / tile_size) * tile_size;
  return aq.parallel_for(sycl::range<1>{a_max}, events, [=](sycl::id<1> id) {
    if (id < array_size)
      VC[id] = static_cast<T2>(VA[id]);
  });
}

template sycl::event copy_n(sycl::queue& aq,
                            const double* restrict VA,
                            size_t array_size,
                            float* restrict VC,
                            const std::vector<sycl::event>& events);

template sycl::event copy_n(sycl::queue& aq,
                            const std::complex<double>* restrict VA,
                            size_t array_size,
                            std::complex<float>* restrict VC,
                            const std::vector<sycl::event>& events);

} // namespace syclBLAS

} // namespace qmcplusplus
