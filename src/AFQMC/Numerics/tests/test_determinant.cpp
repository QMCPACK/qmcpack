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
#include "complex_approx.hpp"
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
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/determinant.hpp"
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

TEST_CASE("determinant_from_getrf", "[Numerics][determinant]")
{
  Alloc<double> alloc{};
  // From scipy.linalg.lapack.dgetrf(x)
  // C order
  Tensor2D<double> x = {{0.8734294027918162, 0.968540662820932, 0.86919454021392},
                        {0.530855691555599, 0.2327283279772907, 0.011398804277429897},
                        {0.4304688182924905, 0.4023513600373648, 0.5226746713525696}};
  // Fortran order
  Tensor2D<double> lu = {{0.968540662820932, 0.2327283279772907, 0.4023513600373648},
                         {0.9017994146450201, 0.32098142161436394, 0.06762859732916715},
                         {0.897426998761559, -0.6151691884001227, 0.20319672719822834}};
  Tensor1D<int> pivot = {1, 1, 2};
  double log_factor   = 0.0;
  double detx         = 0.06317052169675352;
  using ma::determinant_from_getrf;
  double ovlp = determinant_from_getrf(std::get<0>(x.sizes()), lu.origin(), std::get<1>(lu.sizes()), pivot.origin(), log_factor);
  CHECK(ovlp == Approx(detx));
}

TEST_CASE("strided_determinant_from_getrf", "[Numerics][determinant]")
{
  Alloc<double> alloc{};
  // From scipy.linalg.lapack.dgetrf(x)
  // C order
  Tensor2D<double> x = {{0.8734294027918162, 0.968540662820932, 0.86919454021392},
                        {0.530855691555599, 0.2327283279772907, 0.011398804277429897},
                        {0.4304688182924905, 0.4023513600373648, 0.5226746713525696}};
  // Fortran order
  Tensor2D<double> lu = {{0.968540662820932, 0.2327283279772907, 0.4023513600373648},
                         {0.9017994146450201, 0.32098142161436394, 0.06762859732916715},
                         {0.897426998761559, -0.6151691884001227, 0.20319672719822834}};
  Tensor3D<double> lus({3, 3, 3}, 0.0, alloc);
  using std::copy_n;
  for (int i = 0; i < 3; i++)
    copy_n(lu.origin(), lu.num_elements(), lus[i].origin());
  Tensor2D<int> pivot    = {{1, 1, 2}, {1, 1, 2}, {1, 1, 2}};
  Tensor1D<double> ovlps = {0, 0, 0};
  double log_factor      = 0.0;
  double detx            = 0.06317052169675352;
  using ma::strided_determinant_from_getrf;
  strided_determinant_from_getrf(std::get<0>(x.sizes()), lus.origin(), std::get<1>(lu.sizes()), lu.num_elements(), pivot.origin(), std::get<1>(pivot.sizes()),
                                 log_factor, to_address(ovlps.origin()), std::get<0>(lus.sizes()));
  CHECK(ovlps[0] == Approx(detx));
  CHECK(ovlps[1] == Approx(detx));
  CHECK(ovlps[2] == Approx(detx));
}

TEST_CASE("batched_determinant_from_getrf", "[Numerics][determinant]")
{
  Alloc<double> alloc{};
  // From scipy.linalg.lapack.dgetrf(x)
  // C order
  Tensor2D<double> x = {{0.8734294027918162, 0.968540662820932, 0.86919454021392},
                        {0.530855691555599, 0.2327283279772907, 0.011398804277429897},
                        {0.4304688182924905, 0.4023513600373648, 0.5226746713525696}};
  // Fortran order
  Tensor2D<double> lu = {{0.968540662820932, 0.2327283279772907, 0.4023513600373648},
                         {0.9017994146450201, 0.32098142161436394, 0.06762859732916715},
                         {0.897426998761559, -0.6151691884001227, 0.20319672719822834}};
  std::vector<pointer<double>> lu_array;
  using std::copy_n;
  for (int i = 0; i < 3; i++)
    lu_array.emplace_back(lu.origin());
  Tensor2D<int> pivot    = {{1, 1, 2}, {1, 1, 2}, {1, 1, 2}};
  Tensor1D<double> ovlps = {0, 0, 0};
  double log_factor      = 0.0;
  double detx            = 0.06317052169675352;
  using ma::batched_determinant_from_getrf;
  batched_determinant_from_getrf(std::get<0>(x.sizes()), lu_array.data(), std::get<1>(lu.sizes()), pivot.origin(), std::get<1>(pivot.sizes()), log_factor,
                                 to_address(ovlps.origin()), lu_array.size());
  CHECK(ovlps[0] == Approx(detx));
  CHECK(ovlps[1] == Approx(detx));
  CHECK(ovlps[2] == Approx(detx));
}

TEST_CASE("batched_determinant_from_getrf_complex", "[Numerics][determinant]")
{
  Alloc<std::complex<double>> alloc{};
  // From scipy.linalg.lapack.dgetrf(x)
  // C order
  Tensor2D<std::complex<double>> x = {{0.8734294027918162, 0.968540662820932, 0.86919454021392},
                                      {0.530855691555599, 0.2327283279772907, 0.011398804277429897},
                                      {0.4304688182924905, 0.4023513600373648, 0.5226746713525696}};
  // Fortran order
  Tensor2D<std::complex<double>> lu = {{0.968540662820932, 0.2327283279772907, 0.4023513600373648},
                                       {0.9017994146450201, 0.32098142161436394, 0.06762859732916715},
                                       {0.897426998761559, -0.6151691884001227, 0.20319672719822834}};
  std::vector<pointer<std::complex<double>>> lu_array;
  using std::copy_n;
  for (int i = 0; i < 3; i++)
    lu_array.emplace_back(lu.origin());
  Tensor2D<int> pivot                  = {{1, 1, 2}, {1, 1, 2}, {1, 1, 2}};
  Tensor1D<std::complex<double>> ovlps = {0, 0, 0};
  std::complex<double> log_factor      = 0.0;
  std::complex<double> detx            = 0.06317052169675352;
  using ma::batched_determinant_from_getrf;
  batched_determinant_from_getrf(std::get<0>(x.sizes()), lu_array.data(), std::get<1>(lu.sizes()), pivot.origin(), std::get<1>(pivot.sizes()), log_factor,
                                 to_address(ovlps.origin()), lu_array.size());
  CHECK(ovlps[0] == ComplexApprox(detx));
  CHECK(ovlps[1] == ComplexApprox(detx));
  CHECK(ovlps[2] == ComplexApprox(detx));
}

TEST_CASE("determinant", "[Numerics][determinant]")
{
  Alloc<std::complex<double>> alloc{};
  Tensor1D<int> pivot = {0, 0, 0, 0};
  Tensor2D<std::complex<double>> WORK({3, 3}, 0.0, alloc);
  double log_ovlp = 0.0;
  // numpy.random.seed(8)
  // x = numpy.random.random(9).reshape(3,3)
  // z = x + 1j*x
  // det = numpy.linalg.det(z)
  using namespace std::complex_literals;
  Tensor2D<std::complex<double>> matx = {{0.8734294027918162 + 0.8734294027918162i,
                                          0.968540662820932 + 0.968540662820932i, 0.86919454021392 + 0.86919454021392i},
                                         {0.530855691555599 + 0.530855691555599i,
                                          0.2327283279772907 + 0.2327283279772907i,
                                          0.011398804277429897 + 0.011398804277429897i},
                                         {0.4304688182924905 + 0.4304688182924905i,
                                          0.4023513600373648 + 0.4023513600373648i,
                                          0.5226746713525696 + 0.5226746713525696i}};
  std::complex<double> ref(0.12634104339350705, -0.12634104339350705);
  using ma::getrf;
  getrf(matx, pivot, WORK.elements());
  using ma::determinant;
  std::complex<double> lovlp         = 0.0;
  Tensor2D<std::complex<double>> mat = {{0.8734294027918162 + 0.8734294027918162i,
                                         0.968540662820932 + 0.968540662820932i, 0.86919454021392 + 0.86919454021392i},
                                        {0.530855691555599 + 0.530855691555599i,
                                         0.2327283279772907 + 0.2327283279772907i,
                                         0.011398804277429897 + 0.011398804277429897i},
                                        {0.4304688182924905 + 0.4304688182924905i,
                                         0.4023513600373648 + 0.4023513600373648i,
                                         0.5226746713525696 + 0.5226746713525696i}};
  std::complex<double> ovlp          = determinant(mat, pivot, WORK.elements(), lovlp);
  CHECK(std::real(ovlp) == Approx(std::real(ref)));
  CHECK(std::imag(ovlp) == Approx(std::imag(ref)));
}

} // namespace qmcplusplus
