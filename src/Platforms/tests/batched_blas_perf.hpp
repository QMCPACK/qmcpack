//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_BATCHEDBLASPERF_HPP
#define QMCPLUSPLUS_BATCHEDBLASPERF_HPP


#include <algorithm>
#include <chrono>
#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <vector>

#include "AccelBLAS.hpp"
#include "Common/AccelBLASHandle.hpp"
#include "DualAllocatorAliases.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"


namespace qmcplusplus
{
namespace benchmark
{
struct BenchConfig
{
  int warmup_iters = 3; // dry runs
  int iters_hint   = 0; // 0 => auto
  int iters_min    = 10;
  int iters_max    = 50000;
  double target_ms = 200.0; // aim total runtime ~= target_ms
};

struct PerfStats
{
  int iters       = 0;
  double ms_total = 0.0;
  double ms_avg   = 0.0;
};

template<typename Queue, typename Submit>
inline PerfStats run_batched_perf(Queue& queue, Submit&& submit, const BenchConfig& cfg = {})
{
  // Warmup (not timed)
  for (int i = 0; i < cfg.warmup_iters; ++i)
    submit(1);
  queue.sync();

  // Time a given number of repetitions
  auto time_once = [&](int reps) -> double {
    queue.sync();
    auto t0 = std::chrono::high_resolution_clock::now();
    submit(reps);
    queue.sync();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
  };

  int iters = cfg.iters_hint;
  // Determine number of iterations if not given
  if (iters <= 0)
  {
    double ms1 = time_once(100) / 100.0;
    ms1        = std::max(ms1, 0.0001);
    iters      = int(std::ceil(cfg.target_ms / ms1));
    iters      = std::max(iters, cfg.iters_min);
    iters      = std::min(iters, cfg.iters_max);
  }

  // Final timed run
  double ms_total = time_once(iters);
  // Report
  PerfStats st;
  st.iters    = iters;
  st.ms_total = ms_total;
  st.ms_avg   = ms_total / double(iters);

  return st;
}

inline void print_perf_line(const char* op,
                            const char* backend,
                            const char* dtype,
                            const char transa,
                            const char transb,
                            const int M,
                            const int N,
                            const int K,
                            const int batch,
                            const PerfStats& st)
{
  std::cout << "[" << op << "] " << backend << " " << dtype
            << (transa != '\0' ? (std::string(" transa=") + transa) : "")
            << (transb != '\0' ? (std::string(" transb=") + transb) : "") << " M=" << M << " N=" << N
            << (K ? (std::string(" K=") + std::to_string(K)) : "") << " batch=" << batch << " iters=" << std::fixed
            << std::setw(5) << std::setfill(' ') << st.iters << " avg=" << std::fixed << std::setprecision(6)
            << st.ms_avg << " ms"
            << "\n";
}

template<PlatformKind PL>
constexpr const char* backend_name();
template<>
constexpr const char* backend_name<PlatformKind::CUDA>()
{
  return "CUDA";
}
template<>
constexpr const char* backend_name<PlatformKind::SYCL>()
{
  return "SYCL";
}
template<>
constexpr const char* backend_name<PlatformKind::OMPTARGET>()
{
  return "OMPT";
}

template<typename T>
constexpr const char* dtype_name();
template<>
constexpr const char* dtype_name<float>()
{
  return " float";
}
template<>
constexpr const char* dtype_name<double>()
{
  return "double";
}
template<>
constexpr const char* dtype_name<std::complex<float>>()
{
  return " cfloat";
}
template<>
constexpr const char* dtype_name<std::complex<double>>()
{
  return "cdouble";
}

// gemv_batched perf helper
template<PlatformKind PL, typename T>
inline void perf_gemv_batched(const char trans, const int M_b, const int N_b, const int batch_size)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  const int M = (trans == 'T') ? M_b : N_b;
  const int N = (trans == 'T') ? N_b : M_b;

  // per-batch buffers
  std::vector<vec_t> X(batch_size);
  std::vector<vec_t> Y(batch_size);
  std::vector<mat_t> B(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    X[b] = vec_t(N);
    Y[b] = vec_t(M);
    B[b] = mat_t(M_b, N_b);
  }

  // Initialize and transfer to device (outside the timed region)
  for (int b = 0; b < batch_size; ++b)
  {
    for (int i = 0; i < N; ++i)
      X[b][i] = T((i + 1) + (b + 1));
    for (int i = 0; i < M; ++i)
      Y[b][i] = T(0);
    for (int j = 0; j < M_b; ++j)
      for (int i = 0; i < N_b; ++i)
        B[b][j][i] = T((i + 1) * (b + 1) + j);

    X[b].updateTo();
    Y[b].updateTo();
    B[b].updateTo();
  }

  // Pointer arrays for batched API
  Vector<const T*, PinnedDualAllocator<const T*>> Xarr(batch_size);
  Vector<const T*, PinnedDualAllocator<const T*>> Barr(batch_size);
  Vector<T*, PinnedDualAllocator<T*>> Yarr(batch_size);
  Vector<T, PinnedDualAllocator<T>> alpha_arr(batch_size);
  Vector<T, PinnedDualAllocator<T>> beta_arr(batch_size);

  for (int b = 0; b < batch_size; ++b)
  {
    Xarr[b]      = X[b].device_data();
    Barr[b]      = B[b].device_data();
    Yarr[b]      = Y[b].device_data();
    alpha_arr[b] = T(1);
    beta_arr[b]  = T(0);
  }
  Xarr.updateTo();
  Barr.updateTo();
  Yarr.updateTo();
  alpha_arr.updateTo();
  beta_arr.updateTo();

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);

  const int lda  = N_b;
  const int incx = 1;
  const int incy = 1;

  // Lambda for gemv_batched for timing
  auto submit = [&](int reps) {
    for (int r = 0; r < reps; ++r)
    {
      compute::BLAS::gemv_batched(h_blas, trans, N_b, M_b, alpha_arr.device_data(), Barr.device_data(), lda,
                                  Xarr.device_data(), incx, beta_arr.device_data(), Yarr.device_data(), incy,
                                  batch_size);
    }
  };

  // Run the performance test and reports
  PerfStats st = run_batched_perf(queue, submit);
  print_perf_line("gemv_batched", backend_name<PL>(), dtype_name<T>(), trans, /*transb=*/'\0', M_b, N_b, /*K=*/0,
                  batch_size, st);
}


// gemm_batched perf helper
template<PlatformKind PL, typename T>
inline void perf_gemm_batched(const char transa,
                              const char transb,
                              const int M,
                              const int N,
                              const int K,
                              const int batch_size)
{
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Dimensions for row-major storage
  const int a0 = (transa == 'T') ? M : K;
  const int a1 = (transa == 'T') ? K : M;

  const int b0 = (transb == 'T') ? K : N;
  const int b1 = (transb == 'T') ? N : K;

  const int c_rows = N;
  const int c_cols = M;

  // Per-batch buffers
  std::vector<mat_t> A(batch_size), B(batch_size), C(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    A[b] = mat_t(a0, a1);
    B[b] = mat_t(b0, b1);
    C[b] = mat_t(c_rows, c_cols);
  }

  // Initialize and transfer to device (outside the timed region)
  for (int b = 0; b < batch_size; ++b)
  {
    for (int j = 0; j < a0; ++j)
      for (int i = 0; i < a1; ++i)
        A[b][j][i] = T(i * 3 + j * 4);

    for (int j = 0; j < b0; ++j)
      for (int i = 0; i < b1; ++i)
        B[b][j][i] = T(i * 4 + j * 5);

    for (int j = 0; j < c_rows; ++j)
      for (int i = 0; i < c_cols; ++i)
        C[b][j][i] = T(0);

    A[b].updateTo();
    B[b].updateTo();
    C[b].updateTo();
  }

  // Pointer arrays for batched API
  Vector<const T*, PinnedDualAllocator<const T*>> Aarr(batch_size), Barr(batch_size);
  Vector<T*, PinnedDualAllocator<T*>> Carr(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    Aarr[b] = A[b].device_data();
    Barr[b] = B[b].device_data();
    Carr[b] = C[b].device_data();
  }
  Aarr.updateTo();
  Barr.updateTo();
  Carr.updateTo();

  const T alpha(1);
  const T beta(0);

  const int lda = a1;
  const int ldb = b1;
  const int ldc = M;

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);

  // Submit loop for timing
  auto submit = [&](int reps) {
    for (int r = 0; r < reps; ++r)
    {
      compute::BLAS::gemm_batched(h_blas, transa, transb, M, N, K, alpha, Aarr.device_data(), lda, Barr.device_data(),
                                  ldb, beta, Carr.device_data(), ldc, batch_size);
    }
  };

  // Time and print
  PerfStats st = run_batched_perf(queue, submit);
  print_perf_line("gemm_batched", backend_name<PL>(), dtype_name<T>(), transa, transb, M, N, K, batch_size, st);
}

// ger_batched perf helper
template<PlatformKind PL, typename T>
inline void perf_ger_batched(const int M, const int N, const int batch_size)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Per-batch buffers
  std::vector<vec_t> x(batch_size);
  std::vector<vec_t> y(batch_size);
  std::vector<mat_t> A(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    x[b] = vec_t(M);
    y[b] = vec_t(N);
    A[b] = mat_t(M, N);
  }

  // Initialize and transfer to device (outside the timed region)
  for (int b = 0; b < batch_size; ++b)
  {
    for (int i = 0; i < M; ++i)
      x[b][i] = T(i);
    for (int j = 0; j < N; ++j)
      y[b][j] = T(N - j);
    for (int r = 0; r < M; ++r)
      for (int c = 0; c < N; ++c)
        A[b][r][c] = T(0);

    x[b].updateTo();
    y[b].updateTo();
    A[b].updateTo();
  }

  // Pointer arrays for batched API
  Vector<T*, PinnedDualAllocator<T*>> Aarr(batch_size);
  Vector<const T*, PinnedDualAllocator<const T*>> Xarr(batch_size), Yarr(batch_size);
  Vector<T, PinnedDualAllocator<T>> alpha_arr(batch_size);

  for (int b = 0; b < batch_size; ++b)
  {
    Aarr[b]      = A[b].device_data();
    Xarr[b]      = x[b].device_data();
    Yarr[b]      = y[b].device_data();
    alpha_arr[b] = T(1);
  }
  Aarr.updateTo();
  Xarr.updateTo();
  Yarr.updateTo();
  alpha_arr.updateTo();

  const int incx = 1;
  const int incy = 1;
  const int lda  = M;

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);

  // Timed submit loop
  auto submit = [&](int reps) {
    for (int r = 0; r < reps; ++r)
    {
      compute::BLAS::ger_batched(h_blas, M, N, alpha_arr.device_data(), Xarr.device_data(), incx, Yarr.device_data(),
                                 incy, Aarr.device_data(), lda, batch_size);
    }
  };

  PerfStats st = run_batched_perf(queue, submit);
  print_perf_line("ger_batched", backend_name<PL>(), dtype_name<T>(), '\0', '\0', M, N, /*K=*/0, batch_size, st);
}

} // namespace benchmark
} // namespace qmcplusplus

#endif
