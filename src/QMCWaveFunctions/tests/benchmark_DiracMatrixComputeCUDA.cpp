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

#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <catch.hpp>
#include <algorithm>
#include "Configuration.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Fermion/DiracMatrixComputeCUDA.hpp"
#include "QMCWaveFunctions/tests/makeRngSpdMatrix.hpp"
#include "Utilities/for_testing/checkMatrix.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include "Platforms/PinnedAllocator.h"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"

// Legacy CPU inversion for temporary testing
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"


namespace qmcplusplus
{
template<typename T>
using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
template<typename T>
using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
template<typename T>
using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;


TEST_CASE("DiracMatrixComputeCUDA_large_determinants_benchmark_legacy", "[wavefunction][fermion][benchmark]")
{
  int n = 1024;
  int batch_size = 8;
  

  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  DiracMatrixComputeCUDA<double> dmcc(cuda_handles->hstream);

  Matrix<double> mat_spd;
  mat_spd.resize(n, n);
  makeRngSpdMatrix(mat_spd);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<double> mat_a(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a(i, j) = mat_spd(i, j);

  Matrix<double> mat_spd2;
  mat_spd2.resize(n, n);
  makeRngSpdMatrix(mat_spd2);
  // You would hope you could do this
  // OffloadPinnedMatrix<double> mat_a(mat_spd);
  // But you can't
  OffloadPinnedMatrix<double> mat_a2(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_a2(i, j) = mat_spd2(i, j);
  
  OffloadPinnedVector<std::complex<double>> log_values;
  log_values.resize(2);
  OffloadPinnedMatrix<double> inv_mat_a;
  inv_mat_a.resize(n, n);
  OffloadPinnedMatrix<double> inv_mat_a2;
  inv_mat_a2.resize(n, n);

  RefVector<OffloadPinnedMatrix<double>> a_mats{mat_a, mat_a2};
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats{inv_mat_a, inv_mat_a2};

  BENCHMARK_ADVANCED("CUDA 2x1024")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values, {true, true});});
  };
  
  DiracMatrix<double> dmat;
  Matrix<double> inv_mat_test(n, n);
  std::complex<double> det_log_value;
  Matrix<double> inv_mat_test2(n, n);
  std::complex<double> det_log_value2;

  BENCHMARK("LEGACY 2x1024") {
    dmat.invert_transpose(mat_spd, inv_mat_test, det_log_value);
    return dmat.invert_transpose(mat_spd2, inv_mat_test2, det_log_value2);
  }; 
}

TEST_CASE("benchmark_DiracMatrixComputeCUDA_vs_legacy_256_10", "[wavefunction][fermion][benchmark]")
{
  size_t n = 256;
  int batch_size = 10;

  auto cuda_handles = std::make_unique<CUDALinearAlgebraHandles>();
  DiracMatrixComputeCUDA<double> dmcc(cuda_handles->hstream);

  std::vector<Matrix<double>> spd_mats(batch_size,{n,n});
  std::vector<OffloadPinnedMatrix<double>> pinned_spd_mats(batch_size,{n,n});

  for(int im = 0; im < batch_size; ++im)
  {
    makeRngSpdMatrix(spd_mats[im]);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        pinned_spd_mats[im](i, j) = spd_mats[im](i, j);
  }
  
  OffloadPinnedVector<std::complex<double>> log_values(batch_size);
  std::vector<OffloadPinnedMatrix<double>> pinned_inv_mats(batch_size,{n,n});

  auto a_mats = makeRefVector<decltype(pinned_spd_mats)::value_type>(pinned_spd_mats);
  RefVector<OffloadPinnedMatrix<double>> inv_a_mats = makeRefVector<decltype(pinned_inv_mats)::value_type>(pinned_inv_mats);

  std::vector<bool> compute_mask(batch_size, true);
  BENCHMARK_ADVANCED("CUDA 10x256")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {dmcc.mw_invertTranspose(*cuda_handles, a_mats, inv_a_mats, log_values, compute_mask);});
  };

  
  DiracMatrix<double> dmat;
  std::vector<Matrix<double>> inv_mats_test(batch_size,{n,n});;
  std::vector<std::complex<double>> log_values_test(batch_size);

  BENCHMARK_ADVANCED("LEGACY 10x256")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
                    for(int im = 0; im < batch_size; ++im)
                      dmat.invert_transpose(spd_mats[im], inv_mats_test[im], log_values_test[im]);
                  });
  }; 
}

}
