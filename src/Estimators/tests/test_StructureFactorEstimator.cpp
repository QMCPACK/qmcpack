//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "test_StructureFactorEstimator.h"
#include "StructureFactorInput.h"
#include "ValidStructureFactorInput.h"
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
#include "RandomForTest.h"
#include "GenerateRandomParticleSets.h"
#include "Utilities/ProjectData.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include "Utilities/for_testing/checkVector.hpp"
#include <iostream>

namespace qmcplusplus
{

constexpr bool generate_test_data = false;

using Value = QMCTraits::ValueType;

TEST_CASE("StructureFactorEstimator::StructureFactorEstimator", "[estimators]")
{
  using Input = qmcplusplus::testing::ValidStructureFactorInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  ParticleSetPool particle_pool{MinimalParticlePool::make_diamondC_1x1x1(comm)};

  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::getXml(Input::valid::SKALL));
  xmlNodePtr node = doc.getRoot();
  UPtr<StructureFactorInput> sf_in;
  sf_in = std::make_unique<StructureFactorInput>(node);
  StructureFactorEstimator sfe(*sf_in, pset_ions, pset_elec);
  CHECK(sfe.get_name() == "sk1");
}

TEST_CASE("StructureFactorEstimator::Accumulate", "[estimators]")
{
  using Input = qmcplusplus::testing::ValidStructureFactorInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::getXml(Input::valid::SKALL));
  xmlNodePtr node = doc.getRoot();
  UPtr<StructureFactorInput> sf_in;
  sf_in = std::make_unique<StructureFactorInput>(node);

  ParticleSetPool particle_pool{MinimalParticlePool::make_diamondC_1x1x1(comm)};

  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  StructureFactorEstimator sfe(*sf_in, pset_ions, pset_elec);

  std::vector<OperatorEstBase::MCPWalker> walkers;
  int nwalkers = 3;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  std::vector<ParticleSet::ParticlePos> deterministic_rs = {
      {
          {-0.5864843726, 0.7036851048, -0.3468246162},
          {-0.3687131703, -0.0850731656, 1.026543021},
          {1.690787077, 2.543056011, 2.725890636},
          {1.656631947, 1.702485442, 2.588883162},
          {0.741176188, 0.2260864824, 0.9869101048},
          {-2.001264095, 0.5528662801, 1.449588537},
          {1.242099285, 1.414998293, 2.616523981},
          {0.2022707462, 0.9160885215, 1.679067254},
      },
      {
          {0.4825870097, 2.534951687, -1.179359913},
          {0.1491553783, -0.5051442981, -2.215370893},
          {1.080871463, -1.095661044, 2.88787508},
          {4.592846394, 3.424248695, 1.041036963},
          {0.4105701745, 1.359683633, -0.3659555912},
          {0.939719677, -2.190648556, -1.105555773},
          {-0.5726982355, 1.292797208, 2.216959},
          {0.5649974346, 1.815114617, 2.27428937},
      },
      {
          {0.2357321531, -0.3780939579, -0.05520464852},
          {0.247727558, 1.035735369, 1.439867616},
          {0.5015376806, 0.2596193552, 1.314096689},
          {2.564687729, 4.015676022, 1.641964436},
          {-0.7450796366, -1.681605577, 0.8863996863},
          {-0.1126421094, -0.4667637348, 1.119995475},
          {2.070413589, -0.4000400305, 2.256533146},
          {3.036986828, 1.300868988, 1.215788841},
      },
  };
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  std::vector<ParticleSet> psets =
      testing::generateRandomParticleSets(pset_elec, pset_ions, deterministic_rs, nwalkers, generate_test_data);

  auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
  RefVectorWithLeader<OperatorEstBase::MCPWalker> rvwl_walkers(ref_walkers[0], ref_walkers);
  auto ref_psets = makeRefVector<ParticleSet>(psets);
  RefVectorWithLeader<ParticleSet> rvwl_psets(ref_psets[0], ref_psets);

  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);

  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& spomap = wavefunction_pool.getWaveFunction("wavefunction")->getSPOMap();

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  auto ref_wfns = convertUPtrToRefVector(twfcs);

  // These hamiltomians are just pro forma arguments needed to hold off UBSan,
  // StructureFactorEstimator never accesses into them.
  auto hamiltonian_pool  = MinimalHamiltonianPool::makeHamWithEEEI(comm, particle_pool, wavefunction_pool);
  auto& gold_hamiltonian = *(hamiltonian_pool.getPrimary());
  std::vector<UPtr<QMCHamiltonian>> hams(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    hams[iw] = gold_hamiltonian.makeClone(psets[iw], ref_wfns[iw]);

  auto ref_hams = convertUPtrToRefVector(hams);
  RefVectorWithLeader<QMCHamiltonian> rvwl_hams(ref_hams[0], ref_hams);

  auto updateWalker = [](auto& walker, auto& pset_target, auto& trial_wavefunction) {
    pset_target.update(false);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    //pset_target.saveWalker(walker);
  };

  using QMCT = OperatorEstBase::QMCT;
  std::vector<QMCT::RealType> rng_reals(nwalkers * QMCT::DIM * 2);

  if (generate_test_data)
  {
    testing::RandomForTest<QMCT::RealType> rft;
    rft.fillVecRng(rng_reals);
    std::cout << "rng_reals = " << NativePrint(rng_reals) << '\n';
  }
  else
  {
    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
    rng_reals = {
        0.6121701598,  0.120757781,  0.1690697521, 0.3017894626, 0.4360590279, 0.8539064527,
        0.7692624927,  0.2977159917, 0.295325309,  0.3778624535, 0.149162963,  0.2090236992,
        0.02247832529, 0.6495640278, 0.4202244878, 0.6017025113, 0.2386821359, 0.7122730613,
    };
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  auto it_rng_reals = rng_reals.begin();

  for (ParticleSet& pset : psets)
  {
    // gcc will evaluate these iterators in the opposite of the order explicitly imposed if we naively
    // increment the iterators after each value.
    pset.R[0] = ParticleSet::PosType(*it_rng_reals, *(it_rng_reals + 1), *(it_rng_reals + 2));
    it_rng_reals += 3;
    pset.R[1] = ParticleSet::PosType(*it_rng_reals, *(it_rng_reals + 1), *(it_rng_reals + 2));
    it_rng_reals += 3;
  }


  psets[0].mw_update(rvwl_psets, false);
  // for (int iw = 0; iw < nwalkers; ++iw)
  //   updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  CHECK(sfe.getNumKPoints() == 608);


  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  std::remove_cv_t<std::remove_reference_t<decltype(KContainer().getKpts())>> kpoint_lists = {
      {-1, -1, -1}, {-1, 0, 0},   {0, -1, 0},   {0, 0, -1},   {0, 0, 1},    {0, 1, 0},    {1, 0, 0},    {1, 1, 1},
      {-1, -1, 0},  {-1, 0, -1},  {0, -1, -1},  {0, 1, 1},    {1, 0, 1},    {1, 1, 0},    {-2, -1, -1}, {-1, -2, -1},
      {-1, -1, -2}, {-1, 0, 1},   {-1, 1, 0},   {0, -1, 1},   {0, 1, -1},   {1, -1, 0},   {1, 0, -1},   {1, 1, 2},
      {1, 2, 1},    {2, 1, 1},    {-2, -2, -1}, {-2, -1, -2}, {-2, -1, 0},  {-2, 0, -1},  {-1, -2, -2}, {-1, -2, 0},
      {-1, -1, 1},  {-1, 0, -2},  {-1, 1, -1},  {-1, 1, 1},   {0, -2, -1},  {0, -1, -2},  {0, 1, 2},    {0, 2, 1},
      {1, -1, -1},  {1, -1, 1},   {1, 0, 2},    {1, 1, -1},   {1, 2, 0},    {1, 2, 2},    {2, 0, 1},    {2, 1, 0},
      {2, 1, 2},    {2, 2, 1},    {-2, -2, -2}, {-2, 0, 0},   {0, -2, 0},   {0, 0, -2},   {0, 0, 2},    {0, 2, 0},
      {2, 0, 0},    {2, 2, 2},    {-2, -2, 0},  {-2, 0, -2},  {0, -2, -2},  {0, 2, 2},    {2, 0, 2},    {2, 2, 0},
      {-3, -2, -2}, {-3, -1, -1}, {-2, -3, -2}, {-2, -2, -3}, {-2, 0, 1},   {-2, 1, 0},   {-1, -3, -1}, {-1, -1, -3},
      {-1, 0, 2},   {-1, 2, 0},   {0, -2, 1},   {0, -1, 2},   {0, 1, -2},   {0, 2, -1},   {1, -2, 0},   {1, 0, -2},
      {1, 1, 3},    {1, 3, 1},    {2, -1, 0},   {2, 0, -1},   {2, 2, 3},    {2, 3, 2},    {3, 1, 1},    {3, 2, 2},
      {-3, -2, -1}, {-3, -1, -2}, {-2, -3, -1}, {-2, -1, -3}, {-2, -1, 1},  {-2, 1, -1},  {-1, -3, -2}, {-1, -2, -3},
      {-1, -2, 1},  {-1, 1, -2},  {-1, 1, 2},   {-1, 2, 1},   {1, -2, -1},  {1, -1, -2},  {1, -1, 2},   {1, 2, -1},
      {1, 2, 3},    {1, 3, 2},    {2, -1, 1},   {2, 1, -1},   {2, 1, 3},    {2, 3, 1},    {3, 1, 2},    {3, 2, 1},
      {-3, -3, -2}, {-3, -2, -3}, {-3, -1, 0},  {-3, 0, -1},  {-2, -3, -3}, {-2, 1, 1},   {-1, -3, 0},  {-1, -1, 2},
      {-1, 0, -3},  {-1, 2, -1},  {0, -3, -1},  {0, -1, -3},  {0, 1, 3},    {0, 3, 1},    {1, -2, 1},   {1, 0, 3},
      {1, 1, -2},   {1, 3, 0},    {2, -1, -1},  {2, 3, 3},    {3, 0, 1},    {3, 1, 0},    {3, 2, 3},    {3, 3, 2},
      {-3, -3, -3}, {-3, -3, -1}, {-3, -2, 0},  {-3, -1, -3}, {-3, 0, -2},  {-3, 0, 0},   {-2, -3, 0},  {-2, -2, 1},
      {-2, 0, -3},  {-2, 1, -2},  {-1, -3, -3}, {-1, 2, 2},   {0, -3, -2},  {0, -3, 0},   {0, -2, -3},  {0, 0, -3},
      {0, 0, 3},    {0, 2, 3},    {0, 3, 0},    {0, 3, 2},    {1, -2, -2},  {1, 3, 3},    {2, -1, 2},   {2, 0, 3},
      {2, 2, -1},   {2, 3, 0},    {3, 0, 0},    {3, 0, 2},    {3, 1, 3},    {3, 2, 0},    {3, 3, 1},    {3, 3, 3},
      {-4, -2, -2}, {-2, -4, -2}, {-2, -2, -4}, {-2, 0, 2},   {-2, 2, 0},   {0, -2, 2},   {0, 2, -2},   {2, -2, 0},
      {2, 0, -2},   {2, 2, 4},    {2, 4, 2},    {4, 2, 2},    {-4, -3, -2}, {-4, -2, -3}, {-4, -2, -1}, {-4, -1, -2},
      {-3, -4, -2}, {-3, -2, -4}, {-3, -1, 1},  {-3, 1, -1},  {-2, -4, -3}, {-2, -4, -1}, {-2, -3, -4}, {-2, -1, -4},
      {-2, -1, 2},  {-2, 1, 2},   {-2, 2, -1},  {-2, 2, 1},   {-1, -4, -2}, {-1, -3, 1},  {-1, -2, -4}, {-1, -2, 2},
      {-1, 1, -3},  {-1, 1, 3},   {-1, 2, -2},  {-1, 3, 1},   {1, -3, -1},  {1, -2, 2},   {1, -1, -3},  {1, -1, 3},
      {1, 2, -2},   {1, 2, 4},    {1, 3, -1},   {1, 4, 2},    {2, -2, -1},  {2, -2, 1},   {2, -1, -2},  {2, 1, -2},
      {2, 1, 4},    {2, 3, 4},    {2, 4, 1},    {2, 4, 3},    {3, -1, 1},   {3, 1, -1},   {3, 2, 4},    {3, 4, 2},
      {4, 1, 2},    {4, 2, 1},    {4, 2, 3},    {4, 3, 2},    {-4, -3, -3}, {-4, -1, -1}, {-3, -4, -3}, {-3, -3, -4},
      {-3, -3, 0},  {-3, 0, -3},  {-3, 0, 1},   {-3, 1, 0},   {-1, -4, -1}, {-1, -1, -4}, {-1, 0, 3},   {-1, 3, 0},
      {0, -3, -3},  {0, -3, 1},   {0, -1, 3},   {0, 1, -3},   {0, 3, -1},   {0, 3, 3},    {1, -3, 0},   {1, 0, -3},
      {1, 1, 4},    {1, 4, 1},    {3, -1, 0},   {3, 0, -1},   {3, 0, 3},    {3, 3, 0},    {3, 3, 4},    {3, 4, 3},
      {4, 1, 1},    {4, 3, 3},    {-4, -3, -1}, {-4, -1, -3}, {-3, -4, -1}, {-3, -2, 1},  {-3, -1, -4}, {-3, 1, -2},
      {-2, -3, 1},  {-2, 1, -3},  {-1, -4, -3}, {-1, -3, -4}, {-1, 2, 3},   {-1, 3, 2},   {1, -3, -2},  {1, -2, -3},
      {1, 3, 4},    {1, 4, 3},    {2, -1, 3},   {2, 3, -1},   {3, -1, 2},   {3, 1, 4},    {3, 2, -1},   {3, 4, 1},
      {4, 1, 3},    {4, 3, 1},    {-4, -4, -3}, {-4, -3, -4}, {-4, -1, 0},  {-4, 0, -1},  {-3, -4, -4}, {-3, 1, 1},
      {-1, -4, 0},  {-1, -1, 3},  {-1, 0, -4},  {-1, 3, -1},  {0, -4, -1},  {0, -1, -4},  {0, 1, 4},    {0, 4, 1},
      {1, -3, 1},   {1, 0, 4},    {1, 1, -3},   {1, 4, 0},    {3, -1, -1},  {3, 4, 4},    {4, 0, 1},    {4, 1, 0},
      {4, 3, 4},    {4, 4, 3},    {-4, -4, -2}, {-4, -2, -4}, {-4, -2, 0},  {-4, 0, -2},  {-2, -4, -4}, {-2, -4, 0},
      {-2, -2, 2},  {-2, 0, -4},  {-2, 2, -2},  {-2, 2, 2},   {0, -4, -2},  {0, -2, -4},  {0, 2, 4},    {0, 4, 2},
      {2, -2, -2},  {2, -2, 2},   {2, 0, 4},    {2, 2, -2},   {2, 4, 0},    {2, 4, 4},    {4, 0, 2},    {4, 2, 0},
      {4, 2, 4},    {4, 4, 2},    {-4, -4, -4}, {-4, 0, 0},   {0, -4, 0},   {0, 0, -4},   {0, 0, 4},    {0, 4, 0},
      {4, 0, 0},    {4, 4, 4},    {-5, -3, -3}, {-5, -2, -2}, {-4, -4, -1}, {-4, -3, 0},  {-4, -1, -4}, {-4, 0, -3},
      {-3, -5, -3}, {-3, -4, 0},  {-3, -3, -5}, {-3, -3, 1},  {-3, 0, -4},  {-3, 0, 2},   {-3, 1, -3},  {-3, 2, 0},
      {-2, -5, -2}, {-2, -2, -5}, {-2, 0, 3},   {-2, 3, 0},   {-1, -4, -4}, {-1, 3, 3},   {0, -4, -3},  {0, -3, -4},
      {0, -3, 2},   {0, -2, 3},   {0, 2, -3},   {0, 3, -2},   {0, 3, 4},    {0, 4, 3},    {1, -3, -3},  {1, 4, 4},
      {2, -3, 0},   {2, 0, -3},   {2, 2, 5},    {2, 5, 2},    {3, -2, 0},   {3, -1, 3},   {3, 0, -2},   {3, 0, 4},
      {3, 3, -1},   {3, 3, 5},    {3, 4, 0},    {3, 5, 3},    {4, 0, 3},    {4, 1, 4},    {4, 3, 0},    {4, 4, 1},
      {5, 2, 2},    {5, 3, 3},    {-5, -3, -2}, {-5, -2, -3}, {-3, -5, -2}, {-3, -2, -5}, {-3, -1, 2},  {-3, 2, -1},
      {-2, -5, -3}, {-2, -3, -5}, {-2, 1, 3},   {-2, 3, 1},   {-1, -3, 2},  {-1, 2, -3},  {1, -2, 3},   {1, 3, -2},
      {2, -3, -1},  {2, -1, -3},  {2, 3, 5},    {2, 5, 3},    {3, -2, 1},   {3, 1, -2},   {3, 2, 5},    {3, 5, 2},
      {5, 2, 3},    {5, 3, 2},    {-5, -4, -3}, {-5, -3, -4}, {-5, -2, -1}, {-5, -1, -2}, {-4, -5, -3}, {-4, -3, -5},
      {-4, -1, 1},  {-4, 1, -1},  {-3, -5, -4}, {-3, -4, -5}, {-3, 1, 2},   {-3, 2, 1},   {-2, -5, -1}, {-2, -1, -5},
      {-2, -1, 3},  {-2, 3, -1},  {-1, -5, -2}, {-1, -4, 1},  {-1, -2, -5}, {-1, -2, 3},  {-1, 1, -4},  {-1, 1, 4},
      {-1, 3, -2},  {-1, 4, 1},   {1, -4, -1},  {1, -3, 2},   {1, -1, -4},  {1, -1, 4},   {1, 2, -3},   {1, 2, 5},
      {1, 4, -1},   {1, 5, 2},    {2, -3, 1},   {2, 1, -3},   {2, 1, 5},    {2, 5, 1},    {3, -2, -1},  {3, -1, -2},
      {3, 4, 5},    {3, 5, 4},    {4, -1, 1},   {4, 1, -1},   {4, 3, 5},    {4, 5, 3},    {5, 1, 2},    {5, 2, 1},
      {5, 3, 4},    {5, 4, 3},    {-5, -4, -4}, {-5, -4, -2}, {-5, -3, -1}, {-5, -2, -4}, {-5, -1, -3}, {-5, -1, -1},
      {-4, -5, -4}, {-4, -5, -2}, {-4, -4, -5}, {-4, -2, -5}, {-4, -2, 1},  {-4, 0, 1},   {-4, 1, -2},  {-4, 1, 0},
      {-3, -5, -1}, {-3, -2, 2},  {-3, -1, -5}, {-3, 2, -2},  {-2, -5, -4}, {-2, -4, -5}, {-2, -4, 1},  {-2, -3, 2},
      {-2, 1, -4},  {-2, 2, -3},  {-2, 2, 3},   {-2, 3, 2},   {-1, -5, -3}, {-1, -5, -1}, {-1, -3, -5}, {-1, -1, -5},
      {-1, 0, 4},   {-1, 2, 4},   {-1, 4, 0},   {-1, 4, 2},   {0, -4, 1},   {0, -1, 4},   {0, 1, -4},   {0, 4, -1},
      {1, -4, -2},  {1, -4, 0},   {1, -2, -4},  {1, 0, -4},   {1, 1, 5},    {1, 3, 5},    {1, 5, 1},    {1, 5, 3},
      {2, -3, -2},  {2, -2, -3},  {2, -2, 3},   {2, -1, 4},   {2, 3, -2},   {2, 4, -1},   {2, 4, 5},    {2, 5, 4},
      {3, -2, 2},   {3, 1, 5},    {3, 2, -2},   {3, 5, 1},    {4, -1, 0},   {4, -1, 2},   {4, 0, -1},   {4, 2, -1},
      {4, 2, 5},    {4, 4, 5},    {4, 5, 2},    {4, 5, 4},    {5, 1, 1},    {5, 1, 3},    {5, 2, 4},    {5, 3, 1},
      {5, 4, 2},    {5, 4, 4},    {-4, -4, 0},  {-4, 0, -4},  {0, -4, -4},  {0, 4, 4},    {4, 0, 4},    {4, 4, 0},
      {-5, -5, -3}, {-5, -3, -5}, {-5, -2, 0},  {-5, 0, -2},  {-3, -5, -5}, {-3, 2, 2},   {-2, -5, 0},  {-2, -2, 3},
      {-2, 0, -5},  {-2, 3, -2},  {0, -5, -2},  {0, -2, -5},  {0, 2, 5},    {0, 5, 2},    {2, -3, 2},   {2, 0, 5},
      {2, 2, -3},   {2, 5, 0},    {3, -2, -2},  {3, 5, 5},    {5, 0, 2},    {5, 2, 0},    {5, 3, 5},    {5, 5, 3},
      {-5, -5, -4}, {-5, -4, -5}, {-5, -4, -1}, {-5, -1, -4}, {-5, -1, 0},  {-5, 0, -1},  {-4, -5, -5}, {-4, -5, -1},
      {-4, -3, 1},  {-4, -1, -5}, {-4, 1, -3},  {-4, 1, 1},   {-3, -4, 1},  {-3, 1, -4},  {-1, -5, -4}, {-1, -5, 0},
      {-1, -4, -5}, {-1, -1, 4},  {-1, 0, -5},  {-1, 3, 4},   {-1, 4, -1},  {-1, 4, 3},   {0, -5, -1},  {0, -1, -5},
      {0, 1, 5},    {0, 5, 1},    {1, -4, -3},  {1, -4, 1},   {1, -3, -4},  {1, 0, 5},    {1, 1, -4},   {1, 4, 5},
      {1, 5, 0},    {1, 5, 4},    {3, -1, 4},   {3, 4, -1},   {4, -1, -1},  {4, -1, 3},   {4, 1, 5},    {4, 3, -1},
      {4, 5, 1},    {4, 5, 5},    {5, 0, 1},    {5, 1, 0},    {5, 1, 4},    {5, 4, 1},    {5, 4, 5},    {5, 5, 4},
  };
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

  //std::cout << "kpoint_lists = " << NativePrint(sfe.getKLists().kpts) << '\n';

  double tolerance = 0.1;
  {
    INFO("Checking kpoint_lists");
    auto check = checkVector(sfe.getKLists().getKpts(), kpoint_lists, true);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }

  FakeRandom<OHMMS_PRECISION_FULL> rng;

  sfe.accumulate(ref_walkers, ref_psets, ref_wfns, ref_hams, rng);

  testing::StructureFactorAccess sfa;
  auto& sfk_e_e = sfa.getSKElecElec(sfe);
  auto& rhok_e  = sfa.getRhoKElec(sfe);

  if constexpr (generate_test_data)
  {
    std::cout << "sfk_e_e_expected = ";
    std::cout << NativePrint(sfk_e_e) << '\n';
    std::cout << "rhok_e_expected = ";
    std::cout << NativePrint(rhok_e) << '\n';
    //    FAIL_CHECK("Test always fails when generating new test reference data.");
  }

  auto sfk_e_e_expected = StructureFactorTests::getSKElecElec<Value>();
  auto rhok_e_expected  = StructureFactorTests::getRhoKElec<Value>();

  tolerance = 1e-6;
  // With this wider tolerance we can test mixed precision against full precision.
  // 2e-6 is possible if we keep separate sets of mixed test values.
  if constexpr (std::is_same_v<RealAlias<QMCTraits::ValueType>, float>)
    tolerance = 3e-5;
  {
    INFO("In sfk_e_e test");
    auto check = checkVector(sfk_e_e_expected, sfk_e_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL_CHECK(check.result_message); }
  }
  {
    INFO("In rhok_e test");
    auto check = checkVector(rhok_e_expected, rhok_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL_CHECK(check.result_message); }
  }

  PooledData<double> buffer;
  auto bsize = buffer.size();
  std::cout << "Pooled Data before structure factor pack buffer size = " << bsize << '\n';
  sfe.packData(buffer);
  bsize = buffer.size();
  std::cout << "Pooled Data buffer size = " << bsize << '\n';

  buffer.rewind();
  double scale_by = 2.0;
  buffer *= decltype(buffer)::value_type{scale_by};
  sfe.unpackData(buffer);

  auto& sfk_e_e2 = sfa.getSKElecElec(sfe);
  auto& rhok_e2  = sfa.getRhoKElec(sfe);

  Vector<double> sfk_e_e_expected_scaled{sfk_e_e_expected};
  CHECK(sfk_e_e2.size() == sfk_e_e_expected.size());
  sfk_e_e_expected_scaled *= scale_by;
  Vector<std::complex<double>> rhok_e_expected_scaled{rhok_e_expected};
  CHECK(rhok_e2.size() == rhok_e_expected.size());
  rhok_e_expected_scaled *= scale_by;
  {
    INFO("In sfk_e_e pack unpack test");
    auto check = checkVector(sfk_e_e_expected_scaled, sfk_e_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL_CHECK(check.result_message); }
  }
  {
    INFO("In rhok_e pack unpack test");
    auto check = checkVector(rhok_e_expected_scaled, rhok_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL_CHECK(check.result_message); }
  }
}


} // namespace qmcplusplus
