//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** \file
 *  This implements micro benchmarking on DiracMatrixComputeCUDA.hpp
 *  Currently it also benchmarks the same size matrix's and batch sizes
 *  using the Legacy DiracMatrix serially.
 */

#include "catch.hpp"

#include <algorithm>
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeCUDA.hpp"
#include "makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include "Platforms/DualAllocatorAliases.hpp"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"

// Legacy CPU inversion for temporary testing
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"

namespace qmcplusplus
{
template<typename T>
using OffloadPinnedMatrix = Matrix<T, PinnedDualAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, PinnedDualAllocator<T>>;

// Mechanism to pretty print benchmark names.
struct DiracComputeBenchmarkParameters;
std::ostream& operator<<(std::ostream& out, const DiracComputeBenchmarkParameters& dcbmp);
struct DiracComputeBenchmarkParameters
{
  std::string name;
  size_t n;
  int batch_size;
  std::string str()
  {
    std::stringstream stream;
    stream << *this;
    return stream.str();
  }
};

std::ostream& operator<<(std::ostream& out, const DiracComputeBenchmarkParameters& dcbmp)
{
  out << dcbmp.name << " n=" << dcbmp.n << " batch=" << dcbmp.batch_size;
  return out;
}

/** This and other [.benchmark] benchmarks only run if "[benchmark]" is explicitly passed as tag to test.
 */
TEST_CASE("DiracMatrixComputeCUDA_large_determinants_benchmark_legacy_1024_4", "[wavefunction][fermion][.benchmark]")
{
  DiracComputeBenchmarkParameters params;
  params.name       = "Batched CUDA";
  params.n          = 1024;
  params.batch_size = 4;

  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  std::vector<Matrix<double>> spd_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats(params.batch_size, {params.n, params.n});

  qmcplusplus::testing::MakeRngSpdMatrix<double> makeRngSpdMatrix{};
  for (int im = 0; im < params.batch_size; ++im)
  {
    makeRngSpdMatrix(spd_mats[im]);
    for (int i = 0; i < params.n; ++i)
      for (int j = 0; j < params.n; ++j)
        pinned_spd_mats[im](i, j) = spd_mats[im](i, j);
  }

  OffloadPinnedVector<std::complex<double>> log_values(params.batch_size);
  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats(params.batch_size, {params.n, params.n});

  auto a_mats = makeRefVector<const decltype(pinned_spd_mats)::value_type>(pinned_spd_mats);
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats =
      makeRefVector<decltype(pinned_inv_mats)::value_type>(pinned_inv_mats);

  std::vector<bool> compute_mask(params.batch_size, true);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] { dmcc.mw_invertTranspose(cuda_handles, a_mats, inv_a_mats, log_values); });
  };

  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test(params.batch_size, {params.n, params.n});
  ;
  std::vector<std::complex<double>> log_values_test(params.batch_size);

  params.name = "legacy CPU";
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmat.invert_transpose(spd_mats[im], inv_mats_test[im], log_values_test[im]);
    });
  };
}

/** This test will run by default.
 */
TEST_CASE("benchmark_DiracMatrixComputeCUDA_vs_legacy_256_10", "[wavefunction][fermion][benchmark]")
{
  DiracComputeBenchmarkParameters params;
  params.name       = "Batched CUDA";
  params.n          = 256;
  params.batch_size = 10;

  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  std::vector<Matrix<double>> spd_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats(params.batch_size, {params.n, params.n});

  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  for (int im = 0; im < params.batch_size; ++im)
  {
    makeRngSpdMatrix(spd_mats[im]);
    for (int i = 0; i < params.n; ++i)
      for (int j = 0; j < params.n; ++j)
        pinned_spd_mats[im](i, j) = spd_mats[im](i, j);
  }

  OffloadPinnedVector<std::complex<double>> log_values(params.batch_size);
  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats(params.batch_size, {params.n, params.n});

  auto a_mats = makeRefVector<const decltype(pinned_spd_mats)::value_type>(pinned_spd_mats);
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats =
      makeRefVector<decltype(pinned_inv_mats)::value_type>(pinned_inv_mats);

  std::vector<bool> compute_mask(params.batch_size, true);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] { dmcc.mw_invertTranspose(cuda_handles, a_mats, inv_a_mats, log_values); });
  };


  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test(params.batch_size, {params.n, params.n});
  ;
  std::vector<std::complex<double>> log_values_test(params.batch_size);

  params.name = "legacy CPU";
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmat.invert_transpose(spd_mats[im], inv_mats_test[im], log_values_test[im]);
    });
  };
}

/** Only runs if [benchmark] tag is passed.
 */
TEST_CASE("benchmark_DiracMatrixComputeCUDASingle_vs_legacy_256_10", "[wavefunction][fermion][.benchmark]")
{
  DiracComputeBenchmarkParameters params;
  params.name       = "Forced Serial Batched CUDA";
  params.n          = 256;
  params.batch_size = 10;

  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  std::vector<Matrix<double>> spd_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats(params.batch_size, {params.n, params.n});

  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  for (int im = 0; im < params.batch_size; ++im)
  {
    makeRngSpdMatrix(spd_mats[im]);
    for (int i = 0; i < params.n; ++i)
      for (int j = 0; j < params.n; ++j)
        pinned_spd_mats[im](i, j) = spd_mats[im](i, j);
  }

  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedVector<std::complex<double>>> log_values(params.batch_size);
  for (auto& log_value : log_values)
    log_value.resize(1, {0, 0});

  auto a_mats = makeRefVector<decltype(pinned_spd_mats)::value_type>(pinned_spd_mats);
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats =
      makeRefVector<decltype(pinned_inv_mats)::value_type>(pinned_inv_mats);

  std::vector<bool> compute_mask(params.batch_size, true);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmcc.invert_transpose(cuda_handles, pinned_spd_mats[im], pinned_inv_mats[im], log_values[im]);
    });
  };


  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test(params.batch_size, {params.n, params.n});
  ;
  std::vector<std::complex<double>> log_values_test(params.batch_size);

  params.name = "legacy CPU";
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmat.invert_transpose(spd_mats[im], inv_mats_test[im], log_values_test[im]);
    });
  };
}

/** Only runs if [benchmark] tag is passed.
 */
TEST_CASE("benchmark_DiracMatrixComputeCUDASingle_vs_legacy_1024_4", "[wavefunction][fermion][.benchmark]")
{
  DiracComputeBenchmarkParameters params;
  params.name       = "Forced Serial Batched CUDA";
  params.n          = 1024;
  params.batch_size = 4;

  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;

  std::vector<Matrix<double>> spd_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats(params.batch_size, {params.n, params.n});


  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  for (int im = 0; im < params.batch_size; ++im)
  {
    makeRngSpdMatrix(spd_mats[im]);
    for (int i = 0; i < params.n; ++i)
      for (int j = 0; j < params.n; ++j)
        pinned_spd_mats[im](i, j) = spd_mats[im](i, j);
  }

  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats(params.batch_size, {params.n, params.n});
  std::vector<OffloadPinnedVector<std::complex<double>>> log_values(params.batch_size);
  for (auto& log_value : log_values)
    log_value.resize(1, {0, 0});

  auto a_mats = makeRefVector<decltype(pinned_spd_mats)::value_type>(pinned_spd_mats);
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats =
      makeRefVector<decltype(pinned_inv_mats)::value_type>(pinned_inv_mats);

  std::vector<bool> compute_mask(params.batch_size, true);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmcc.invert_transpose(cuda_handles, pinned_spd_mats[im], pinned_inv_mats[im], log_values[im]);
    });
  };


  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test(params.batch_size, {params.n, params.n});
  ;
  std::vector<std::complex<double>> log_values_test(params.batch_size);

  params.name = "legacy CPU";

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int im = 0; im < params.batch_size; ++im)
        dmat.invert_transpose(spd_mats[im], inv_mats_test[im], log_values_test[im]);
    });
  };
}


} // namespace qmcplusplus
