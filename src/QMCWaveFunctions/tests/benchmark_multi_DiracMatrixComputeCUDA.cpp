//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** \file
 *  This implements micro benchmarking on DiracMatrixComputeCUDA.hpp
 *  showing the importance of multiple threads/streams.
 */

#include "catch.hpp"

#include <algorithm>
#include "Configuration.h"
#include "Concurrency/ParallelExecutor.hpp"
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

class PerThreadResources
{
public:
  PerThreadResources(const DiracComputeBenchmarkParameters& params)
      : spd_mats(params.batch_size, {params.n, params.n}),
        pinned_spd_mats(params.batch_size, {params.n, params.n}),
        log_values(params.batch_size),
        pinned_inv_mats(params.batch_size, {params.n, params.n})
  {
    std::cout << "cuda stream: " << cuda_handles.getStream() << '\n';
  }
  CUDALinearAlgebraHandles cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;
  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  std::vector<Matrix<double>> spd_mats;
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats;
  OffloadPinnedVector<std::complex<double>> log_values;
  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats;
  RefVector<const OffloadPinnedMatrix<double>> a_mats;
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats;
};

class SingleStreamPerThreadResources
{
public:
  SingleStreamPerThreadResources(const DiracComputeBenchmarkParameters& params,
                                 CUDALinearAlgebraHandles& cuda_handles_single)
      : cuda_handles(cuda_handles_single),
        spd_mats(params.batch_size, {params.n, params.n}),
        pinned_spd_mats(params.batch_size, {params.n, params.n}),
        log_values(params.batch_size),
        pinned_inv_mats(params.batch_size, {params.n, params.n})
  {}
  CUDALinearAlgebraHandles& cuda_handles;
  DiracMatrixComputeCUDA<double> dmcc;
  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  std::vector<Matrix<double>> spd_mats;
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats;
  OffloadPinnedVector<std::complex<double>> log_values;
  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats;
  RefVector<const OffloadPinnedMatrix<double>> a_mats;
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats;
};

class LegacyResources
{
public:
  LegacyResources(const DiracComputeBenchmarkParameters& params)
      : spd_mats(params.batch_size, {params.n, params.n}),
        inv_mats_test(params.batch_size, {params.n, params.n}),
        log_values_test(params.batch_size)
  {}
  testing::MakeRngSpdMatrix<double> makeRngSpdMatrix;
  std::vector<Matrix<double>> spd_mats;
  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test;
  std::vector<std::complex<double>> log_values_test;
};

/** This test will run by default.
 */

TEST_CASE("benchmark_multi_DiracMatrixComputeCUDA_vs_legacy_256_8", "[wavefunction][fermion][benchmark]")
{
  DiracComputeBenchmarkParameters params;
  params.name       = "Batched CUDA";
  params.n          = 256;
  params.batch_size = 4;

  int n_threads = 4;


  auto setup = [&params](int task_id, std::vector<PerThreadResources>& resources) {
    auto& res = resources[task_id];
    for (int im = 0; im < params.batch_size; ++im)
    {
      res.makeRngSpdMatrix(res.spd_mats[im]);
      for (int i = 0; i < params.n; ++i)
        for (int j = 0; j < params.n; ++j)
          res.pinned_spd_mats[im](i, j) = res.spd_mats[im](i, j);

      if (!(res.a_mats.size() > 0))
      {
        res.a_mats     = makeRefVector<const decltype(res.pinned_spd_mats)::value_type>(res.pinned_spd_mats);
        res.inv_a_mats = makeRefVector<decltype(res.pinned_inv_mats)::value_type>(res.pinned_inv_mats);
      }
    }
  };

  auto invert = [](int task_id, std::vector<PerThreadResources>& resources) {
    auto& res = resources[task_id];
    res.dmcc.mw_invertTranspose(res.cuda_handles, res.a_mats, res.inv_a_mats, res.log_values);
  };

  ParallelExecutor<> task_pool;
  std::vector<PerThreadResources> resources;
  for (int ir = 0; ir < n_threads; ++ir)
    resources.emplace_back(params);


  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    task_pool(n_threads, setup, resources);
    meter.measure([&] { task_pool(n_threads, invert, resources); });
  };

  CUDALinearAlgebraHandles cuda_handles;
  std::vector<SingleStreamPerThreadResources> single_stream_resources;
  for (int ir = 0; ir < n_threads; ++ir)
    single_stream_resources.emplace_back(params, cuda_handles);
  params.name              = "single stream";
  auto setup_single_stream = [&params](int task_id, std::vector<SingleStreamPerThreadResources>& resources) {
    auto& res = resources[task_id];
    for (int im = 0; im < params.batch_size; ++im)
    {
      res.makeRngSpdMatrix(res.spd_mats[im]);
      for (int i = 0; i < params.n; ++i)
        for (int j = 0; j < params.n; ++j)
          res.pinned_spd_mats[im](i, j) = res.spd_mats[im](i, j);

      res.a_mats     = makeRefVector<const decltype(res.pinned_spd_mats)::value_type>(res.pinned_spd_mats);
      res.inv_a_mats = makeRefVector<decltype(res.pinned_inv_mats)::value_type>(res.pinned_inv_mats);
    }
  };

  auto invert_single_stream = [](int task_id, std::vector<SingleStreamPerThreadResources>& resources) {
    auto& res = resources[task_id];
    res.dmcc.mw_invertTranspose(res.cuda_handles, res.a_mats, res.inv_a_mats, res.log_values);
  };

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    task_pool(n_threads, setup_single_stream, single_stream_resources);
    meter.measure([&] { task_pool(n_threads, invert_single_stream, single_stream_resources); });
  };

  params.name = "single thread";

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    task_pool(n_threads, setup_single_stream, single_stream_resources);
    meter.measure([&] {
      for (int ir = 0; ir < n_threads; ++ir)
        invert_single_stream(ir, single_stream_resources);
    });
  };

  params.name = "legacy CPU";

  std::vector<LegacyResources> legacy_resources;
  for (int ir = 0; ir < n_threads; ++ir)
    legacy_resources.emplace_back(params);

  auto setup_legacy = [&params](int task_id, std::vector<LegacyResources>& resources) {
    auto& res = resources[task_id];
    for (int im = 0; im < params.batch_size; ++im)
    {
      res.makeRngSpdMatrix(res.spd_mats[im]);
    }
  };

  auto invert_legacy = [params](int task_id, std::vector<LegacyResources>& resources) {
    auto& res = resources[task_id];
    for (int im = 0; im < params.batch_size; ++im)
      res.dmat.invert_transpose(res.spd_mats[im], res.inv_mats_test[im], res.log_values_test[im]);
  };

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    task_pool(n_threads, setup_legacy, legacy_resources);
    meter.measure([&] { task_pool(n_threads, invert_legacy, legacy_resources); });
  };
}


} // namespace qmcplusplus
