//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//
// File created by: Youngjun Lee, leey@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <DualAllocatorAliases.hpp>
#include <AccelBLAS.hpp>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{

template<PlatformKind PL>
constexpr std::string_view backend_name()
{
  if constexpr (PL == PlatformKind::CUDA)
    return "CUDA";
  else if constexpr (PL == PlatformKind::SYCL)
    return "SYCL";
  else if constexpr (PL == PlatformKind::OMPTARGET)
    return "OMPT";
}

template<typename T>
constexpr std::string_view dtype_name()
{
  if constexpr (std::is_same_v<T, float>)
    return "f32";
  else if constexpr (std::is_same_v<T, double>)
    return "f64";
  else if constexpr (std::is_same_v<T, std::complex<float>>)
    return "c32";
  else if constexpr (std::is_same_v<T, std::complex<double>>)
    return "c64";
}

template<PlatformKind PL, typename T>
std::string make_bench_name(std::string_view op,
                            const int M_b,
                            const int N_b,
                            const int K_b,
                            const int batch_size,
                            const char transa,
                            const char transb)
{
  std::stringstream ss;
  ss << '[' << backend_name<PL>() << '/' << dtype_name<T>() << ']';
  ss << ' ' << op << '_' << M_b << 'x' << N_b;
  if (K_b)
    ss << 'x' << K_b;
  if (transa)
    ss << '_' << transa;
  if (transb)
    ss << transb;
  ss << " bs" << batch_size;
  return ss.str();
}

// gemm_batched benchmark helper
template<PlatformKind PL, typename T>
inline void bench_gemm_batched(const char transa,
                              const char transb,
                              const int M_b,
                              const int N_b,
                              const int K_b,
                              const int batch_size)
{
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Dimensions for row-major storage
  const int a0 = (transa == 'T') ? M_b : K_b;
  const int a1 = (transa == 'T') ? K_b : M_b;

  const int b0 = (transb == 'T') ? K_b : N_b;
  const int b1 = (transb == 'T') ? N_b : K_b;

  const int c_rows = N_b;
  const int c_cols = M_b;

  // Per-batch buffers
  std::vector<mat_t> A(batch_size);
  std::vector<mat_t> B(batch_size);
  std::vector<mat_t> C(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    A[b] = mat_t(a0, a1);
    B[b] = mat_t(b0, b1);
    C[b] = mat_t(c_rows, c_cols);
  }

  // Initialize and transfer to device
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
  Vector<const T*, PinnedDualAllocator<const T*>> Aarr(batch_size);
  Vector<const T*, PinnedDualAllocator<const T*>> Barr(batch_size);
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
  const int ldc = M_b;

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);

  const std::string bench_name = make_bench_name<PL, T>("gemm_batched", M_b, N_b, K_b, batch_size, transa, transb);

  const int iters = 50;
  BENCHMARK_ADVANCED(bench_name.c_str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int i = 0; i < iters; ++i)
      {
        compute::BLAS::gemm_batched(h_blas, transa, transb, M_b, N_b, K_b, alpha, Aarr.device_data(), lda,
                                    Barr.device_data(), ldb, beta, Carr.device_data(), ldc, batch_size);
        queue.sync();
      }
    });
  };
}


// gemv_batched benchmark helper
template<PlatformKind PL, typename T>
inline void bench_gemv_batched(const char trans, const int M_b, const int N_b, const int batch_size)
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

  // Initialize and transfer to device
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

  const std::string bench_name =
      make_bench_name<PL, T>("gemv_batched", M_b, N_b, /*K_b=*/0, batch_size, trans, /*transb=*/'\0');

  const int iters = 100;
  BENCHMARK_ADVANCED(bench_name.c_str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int i = 0; i < iters; ++i)
      {
        compute::BLAS::gemv_batched(h_blas, trans, N_b, M_b, alpha_arr.device_data(), Barr.device_data(), lda,
                                    Xarr.device_data(), incx, beta_arr.device_data(), Yarr.device_data(), incy,
                                    batch_size);
        queue.sync();
      }
    });
  };
}


// ger_batched benchmark helper
template<PlatformKind PL, typename T>
inline void bench_ger_batched(const int M_b, const int N_b, const int batch_size)
{
  using vec_t = Vector<T, PinnedDualAllocator<T>>;
  using mat_t = Matrix<T, PinnedDualAllocator<T>>;

  // Per-batch buffers
  std::vector<vec_t> x(batch_size);
  std::vector<vec_t> y(batch_size);
  std::vector<mat_t> A(batch_size);
  for (int b = 0; b < batch_size; ++b)
  {
    x[b] = vec_t(M_b);
    y[b] = vec_t(N_b);
    A[b] = mat_t(M_b, N_b);
  }

  // Initialize and transfer to device
  for (int b = 0; b < batch_size; ++b)
  {
    for (int i = 0; i < M_b; ++i)
      x[b][i] = T(i);
    for (int j = 0; j < N_b; ++j)
      y[b][j] = T(N_b - j);
    for (int r = 0; r < M_b; ++r)
      for (int c = 0; c < N_b; ++c)
        A[b][r][c] = T(0);

    x[b].updateTo();
    y[b].updateTo();
    A[b].updateTo();
  }

  // Pointer arrays for batched API
  Vector<T*, PinnedDualAllocator<T*>> Aarr(batch_size);
  Vector<const T*, PinnedDualAllocator<const T*>> Xarr(batch_size);
  Vector<const T*, PinnedDualAllocator<const T*>> Yarr(batch_size);
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
  const int lda  = M_b;

  compute::Queue<PL> queue;
  compute::BLASHandle<PL> h_blas(queue);

  const std::string bench_name =
      make_bench_name<PL, T>("ger_batched", M_b, N_b, /*K_b=*/0, batch_size, /*transa=*/'\0', /*transb=*/'\0');

  const int iters = 100;
  BENCHMARK_ADVANCED(bench_name.c_str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int i = 0; i < iters; ++i)
      {
        compute::BLAS::ger_batched(h_blas, M_b, N_b, alpha_arr.device_data(), Xarr.device_data(), incx,
                                   Yarr.device_data(), incy, Aarr.device_data(), lda, batch_size);
        queue.sync();
      }
    });
  };
}


template<PlatformKind PL>
void benchmark_gemm_cases()
{
  const int M_b   = 29;
  const int N_b   = 31;
  const int K_b   = 23;
  const int BATCH = 64;

  // Batched Test
  bench_gemm_batched<PL, float>('N', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, double>('N', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, float>('T', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, double>('T', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, float>('N', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, double>('N', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, float>('T', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, double>('T', 'T', M_b, N_b, K_b, BATCH);
#if defined(QMC_COMPLEX)
  bench_gemm_batched<PL, std::complex<float>>('N', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<double>>('N', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<float>>('T', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<double>>('T', 'N', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<float>>('N', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<double>>('N', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<float>>('T', 'T', M_b, N_b, K_b, BATCH);
  bench_gemm_batched<PL, std::complex<double>>('T', 'T', M_b, N_b, K_b, BATCH);
#endif
}

template<PlatformKind PL>
void benchmark_gemv_cases()
{
  const int M_b   = 137;
  const int N_b   = 79;
  const int BATCH = 64;

  // Batched Test
  bench_gemv_batched<PL, float>('N', M_b, N_b, BATCH);
  bench_gemv_batched<PL, double>('N', M_b, N_b, BATCH);
  bench_gemv_batched<PL, float>('T', M_b, N_b, BATCH);
  bench_gemv_batched<PL, double>('T', M_b, N_b, BATCH);
#if defined(QMC_COMPLEX)
  bench_gemv_batched<PL, std::complex<float>>('N', M_b, N_b, BATCH);
  bench_gemv_batched<PL, std::complex<double>>('N', M_b, N_b, BATCH);
  bench_gemv_batched<PL, std::complex<float>>('T', M_b, N_b, BATCH);
  bench_gemv_batched<PL, std::complex<double>>('T', M_b, N_b, BATCH);
#endif
}

template<PlatformKind PL>
void benchmark_ger_cases()
{
  const int M_b   = 137;
  const int N_b   = 79;
  const int BATCH = 64;

  // Batched Test
  bench_ger_batched<PL, float>(M_b, N_b, BATCH);
  bench_ger_batched<PL, double>(M_b, N_b, BATCH);
#if defined(QMC_COMPLEX)
  bench_ger_batched<PL, std::complex<float>>(M_b, N_b, BATCH);
  bench_ger_batched<PL, std::complex<double>>(M_b, N_b, BATCH);
#endif
}


TEST_CASE("AccelBLAS batched benchmark", "[BLAS][benchmark]")
{
#if defined(ENABLE_CUDA)
  SECTION("gemm_batched<PlatformKind::CUDA>") { benchmark_gemm_cases<PlatformKind::CUDA>(); }
#endif
#if defined(ENABLE_SYCL)
  SECTION("gemm_batched<PlatformKind::SYCL>") { benchmark_gemm_cases<PlatformKind::SYCL>(); }
#endif
#if defined(ENABLE_OFFLOAD)
  SECTION("gemm_batched<PlatformKind::OMPTARGET>") { benchmark_gemm_cases<PlatformKind::OMPTARGET>(); }
#endif

#if defined(ENABLE_CUDA)
  SECTION("gemv_batched<PlatformKind::CUDA>") { benchmark_gemv_cases<PlatformKind::CUDA>(); }
#endif
#if defined(ENABLE_SYCL)
  SECTION("gemv_batched<PlatformKind::SYCL>") { benchmark_gemv_cases<PlatformKind::SYCL>(); }
#endif
#if defined(ENABLE_OFFLOAD)
  SECTION("gemv_batched<PlatformKind::OMPTARGET>") { benchmark_gemv_cases<PlatformKind::OMPTARGET>(); }
#endif

#if defined(ENABLE_CUDA)
  SECTION("ger_batched<PlatformKind::CUDA>") { benchmark_ger_cases<PlatformKind::CUDA>(); }
#endif
#if defined(ENABLE_SYCL)
  SECTION("ger_batched<PlatformKind::SYCL>") { benchmark_ger_cases<PlatformKind::SYCL>(); }
#endif
#if defined(ENABLE_OFFLOAD)
  SECTION("ger_batched<PlatformKind::OMPTARGET>") { benchmark_ger_cases<PlatformKind::OMPTARGET>(); }
#endif
}

} // namespace qmcplusplus
