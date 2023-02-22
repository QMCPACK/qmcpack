//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Fionn D. Malone, malone14@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    exit(0);         \
  }

#include <vector>


#include "AFQMC/config.h"
#include "AFQMC/config.0.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Numerics/batched_operations.hpp"
#include "AFQMC/Matrix/tests/matrix_helpers.h"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"


using boost::multi::array;
using boost::multi::array_ref;
using boost::multi::iextensions;
using std::copy_n;

namespace qmcplusplus
{
using namespace afqmc;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
template<typename T>
using Alloc = device::device_allocator<T>;
#else
template<typename T>
using Alloc = std::allocator<T>;
#endif
template<typename T>
using pointer = typename Alloc<T>::pointer;


template<typename T>
using Tensor1D = array<T, 1, Alloc<T>>;
template<typename T>
using Tensor2D = array<T, 2, Alloc<T>>;
template<typename T>
using Tensor3D = array<T, 3, Alloc<T>>;
template<typename T>
using Tensor4D = array<T, 4, Alloc<T>>;
template<typename T>
using Tensor5D = array<T, 5, Alloc<T>>;

template<typename T>
using Tensor2D_ref = array_ref<T, 2, pointer<T>>;
template<typename T>
using Tensor1D_ref = array_ref<T, 1, pointer<T>>;

template<typename T>
void create_data(std::vector<T>& buffer, T scale)
{
  T count = T(0);
  for (int i = 0; i < buffer.size(); i++)
  {
    buffer[i] = count;
    count += T(1) / scale;
  }
}

TEST_CASE("adotpby", "[Numerics][ma_blas_extensions]")
{
  size_t n = 1025;
  SECTION("double")
  {
    Alloc<double> alloc{};
    Tensor1D<double> y = {0.0};
    Tensor1D<double> a(n);
    Tensor1D<double> b(n);
    double alpha = 0.5;
    double beta  = 0.0;
    for (int i = 0; i < a.size(); i++)
    {
      a[i] = 0.1;
      b[i] = 0.1;
    }
    using ma::adotpby;
    adotpby(alpha, a, b, beta, y.origin());
    CHECK(y[0] == Approx(0.5 * a.size() * 0.01));
  }
  SECTION("complex")
  {
    Alloc<std::complex<double>> alloc{};
    Tensor1D<std::complex<double>> y = {0.0};
    Tensor1D<std::complex<double>> a(n);
    Tensor1D<std::complex<double>> b(n);
    std::complex<double> alpha = 0.5;
    std::complex<double> beta  = 0.0;
    std::vector<std::complex<double>> y_cpu(1);
    for (int i = 0; i < a.size(); i++)
    {
      a[i] = ComplexType(0.1, 0.1);
      b[i] = ComplexType(0.1, -0.1);
    }
    using ma::adotpby;
    adotpby(alpha, a, b, beta, y.origin());
    copy_n(y.data(), y.size(), y_cpu.data());
    CHECK(std::real(y_cpu[0]) == Approx(a.size() * 0.01));
    CHECK(std::imag(y_cpu[0]) == Approx(0.00));
  }
}

TEST_CASE("axty", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  Tensor1D<double> y   = {1.0, 2.0, 3.0, 4.0};
  Tensor1D<double> x   = {0.1, 5.2, 88.4, 0.001};
  Tensor1D<double> ref = {3.33000e-01, 3.46320e+01, 8.83116e+02, 1.33200e-02};
  double alpha         = 3.33;
  using ma::axty;
  axty(alpha, x, y);
  verify_approx(y, ref);
}

TEST_CASE("axty2D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  Tensor2D<double> y   = {{1.0, 2.0, 3.0, 4.0}, {1.0, 2.0, 3.0, 4.0}};
  Tensor2D<double> x   = {{0.1, 5.2, 88.4, 0.001}, {0.1, 5.2, 88.4, 0.001}};
  Tensor2D<double> ref = {{3.33000e-01, 3.46320e+01, 8.83116e+02, 1.33200e-02},
                          {3.33000e-01, 3.46320e+01, 8.83116e+02, 1.33200e-02}};
  double alpha         = 3.33;
  using ma::axty;
  axty(alpha, x, y);
  verify_approx(y, ref);
}

TEST_CASE("acAxpbB", "[Numerics][ma_blas_extensions]")
{
  Alloc<ComplexType> alloc{};
  Tensor2D<ComplexType> z = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
  ComplexType val         = ComplexType(2.33, -0.24);
  Tensor2D<ComplexType> y = {{val, val, val}, {val, val, val}, {val, val, val}};
  Tensor1D<ComplexType> x = {-val, val, -val};
  ComplexType alpha       = 3.33;
  ComplexType beta        = 0.33;
  using ma::acAxpbB;
  acAxpbB(alpha, y, x, beta, z);
  ComplexType r1            = -17.940045;
  ComplexType r2            = 18.600045;
  Tensor2D<ComplexType> ref = {{r1, r2, r1}, {r1, r2, r1}, {r1, r2, r1}};
  verify_approx(z, ref);
}

TEST_CASE("adiagApy", "[Numerics][ma_blas_extensions]")
{
  Alloc<ComplexType> alloc{};
  ComplexType val         = ComplexType(2.33, -0.24);
  Tensor2D<ComplexType> A = {{val, val, val}, {val, val, val}, {val, val, val}};
  Tensor1D<ComplexType> y = {-val, val, -val};
  ComplexType alpha       = 3.33;
  using ma::adiagApy;
  adiagApy(alpha, A, y);
  ComplexType r1            = ComplexType(5.4289, -0.5592);
  ComplexType r2            = ComplexType(10.0889, -1.0392);
  Tensor1D<ComplexType> ref = {r1, r2, r1};
  verify_approx(y, ref);
}

TEST_CASE("sum", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  Tensor1D<double> y = {1, 2, 3};
  using ma::sum;
  double res = sum(y);
  CHECK(res == Approx(6.0));
}

TEST_CASE("sum2D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(3 * 3);
  Tensor2D<double> y({3, 3}, alloc);
  create_data(buffer, 1.0);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::sum;
  double res = sum(y);
  CHECK(res == Approx(36.0));
}

TEST_CASE("sum3D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(3 * 3 * 3);
  Tensor3D<double> y({3, 3, 3}, alloc);
  create_data(buffer, 1.0);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::sum;
  double res = sum(y);
  CHECK(res == Approx(351.0));
}

TEST_CASE("sum4D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(3 * 3 * 3 * 3);
  Tensor4D<double> y({3, 3, 3, 3}, alloc);
  create_data(buffer, 1.0);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::sum;
  double res = sum(y);
  CHECK(res == Approx(3240.0));
}

TEST_CASE("zero_complex_part", "[Numerics][ma_blas_extensions]")
{
  Alloc<ComplexType> alloc{};
  ComplexType val           = ComplexType(1.0, -1.0);
  Tensor1D<ComplexType> y   = {val, val, val};
  Tensor1D<ComplexType> res = {1.0, 1.0, 1.0};
  using ma::zero_complex_part;
  zero_complex_part(y);
  verify_approx(y, res);
}

TEST_CASE("set_identity2D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(3 * 3);
  Tensor2D<double> y({3, 3}, alloc);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::set_identity;
  set_identity(y);
  CHECK(y[0][0] == Approx(1.0));
  CHECK(y[1][1] == Approx(1.0));
  CHECK(y[2][2] == Approx(1.0));
}

TEST_CASE("set_identity3D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(3 * 3 * 3);
  Tensor3D<double> y({3, 3, 3}, alloc);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::set_identity;
  set_identity(y);
  CHECK(y[0][0][0] == Approx(1.0));
  CHECK(y[1][1][1] == Approx(1.0));
  CHECK(y[2][2][2] == Approx(1.0));
}

TEST_CASE("fill2D", "[Numerics][ma_blas_extensions]")
{
  Alloc<double> alloc{};
  std::vector<double> buffer(2 * 2);
  Tensor2D<double> y({2, 2}, alloc);
  copy_n(buffer.data(), buffer.size(), y.origin());
  using ma::fill;
  fill(y, 2.0);
  Tensor2D<double> ref = {{2.0, 2.0}, {2.0, 2.0}};
  verify_approx(y, ref);
}

TEST_CASE("get_diagonal_strided", "[Numerics][ma_blas_extensions]")
{
  Alloc<ComplexType> alloc{};
  int nk = 2;
  int ni = 3;
  int nj = 3;
  Tensor2D<ComplexType> A({nk, ni}, 0.0, alloc);
  Tensor3D<ComplexType> B({nk, ni, nj}, ComplexType(1.0, -3.0), alloc);
  B[0][0][0] = ComplexType(1.0);
  B[0][2][2] = ComplexType(0, -1.0);
  B[1][0][0] = ComplexType(1.0);
  B[1][2][2] = ComplexType(0, -1.0);
  using ma::get_diagonal_strided;
  get_diagonal_strided(B, A);
  Tensor2D<ComplexType> ref = {{ComplexType(1.0), ComplexType(1.0, -3.0), ComplexType(0, -1.0)},
                               {ComplexType(1.0), ComplexType(1.0, -3.0), ComplexType(0, -1.0)}};
  verify_approx(A, ref);
}

} // namespace qmcplusplus
