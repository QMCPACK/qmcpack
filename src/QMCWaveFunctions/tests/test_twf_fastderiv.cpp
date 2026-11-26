//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

/** @file test_twf_fastderiv.cpp
 *  @brief Main unit tests of the TWFFastDerivWrapper API.
 *
 *  Tests here call TWFFastDerivWrapper member functions directly on small hand-built
 *  matrices and check their behavior.
 *
 *  This is in contrast to test_twf_fastderiv_raii.cpp, which only tests how the wrapper's
 *  resources are created, acquired and released when it is managed by TrialWaveFunction.
 */

#include <catch2/catch_test_macros.hpp>

#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "QMCWaveFunctions/tests/FakeSPO.h"

namespace qmcplusplus
{
// With no virtual orbitals (norb == nelec) buildIntermediates, computeMDDerivatives_Obs,
// and computeMDDerivatives_dmu must handle nvirt == 0 without issuing invalid GEMM calls.
TEST_CASE("TWFFastDerivWrapper buildIntermediates without virtual orbitals", "[wavefunction]")
{
  using IndexType   = TWFFastDerivWrapper::IndexType;
  using ValueMatrix = TWFFastDerivWrapper::ValueMatrix;
  using ValueType   = TWFFastDerivWrapper::ValueType;
  using ValueVector = TWFFastDerivWrapper::ValueVector;

  const int n = 2;
  std::vector<ValueMatrix> Minv(1, ValueMatrix(n, n)), B(1, ValueMatrix(n, n)), M(1, ValueMatrix(n, n)),
      X(1, ValueMatrix(n, n)), Minv_B(1, ValueMatrix(n, n)), Minv_Mv(1, ValueMatrix(n, 0));
  Minv[0] = B[0] = M[0] = 1.0;

  auto spo = std::make_unique<FakeSPO<ValueType>>(n);
  SimulationCell simulation_cell;
  ParticleSet P(simulation_cell);

  TWFFastDerivWrapper twf;
  twf.addGroup(P, 0, spo.get());
  twf.buildIntermediates(Minv, B, M, X, Minv_B, Minv_Mv);

  MultiDiracDeterminant mdd(std::move(spo), false, 0, n);
  std::vector<size_t> occ(n, 1);
  std::vector<ci_configuration2> configs{ci_configuration2(occ)};
  std::vector<size_t> c2nodes_unsorted{0};
  std::vector<size_t> c2nodes_sorted(1);
  mdd.createDetData(0, configs, c2nodes_unsorted, c2nodes_sorted);

  const std::vector<IndexType> mdd_spo_ids{0};
  const std::vector<const MultiDiracDeterminant*> mdds{&mdd};

  std::vector<ValueVector> dvals_O(1, ValueVector(1));
  dvals_O[0][0] = 1.0;
  twf.computeMDDerivatives_Obs(Minv_Mv, Minv_B, mdd_spo_ids, mdds, dvals_O);
  CHECK(dvals_O[0][0] == ValueType(0.0));

  std::vector<ValueMatrix> Minv_dM(1, ValueMatrix(n, n)), Minv_dB(1, ValueMatrix(n, n));
  Minv_dM[0] = Minv_dB[0] = 1.0;
  std::vector<ValueVector> dvals_dmu_O(1, ValueVector(1)), dvals_dmu(1, ValueVector(1));
  dvals_dmu_O[0][0] = 1.0;
  dvals_dmu[0][0]   = 1.0;
  twf.computeMDDerivatives_dmu(Minv_Mv, Minv_B, Minv_dM, Minv_dB, mdd_spo_ids, mdds, dvals_dmu_O, dvals_dmu);
  CHECK(dvals_dmu_O[0][0] == ValueType(0.0));
  CHECK(dvals_dmu[0][0] == ValueType(0.0));
}
} // namespace qmcplusplus
