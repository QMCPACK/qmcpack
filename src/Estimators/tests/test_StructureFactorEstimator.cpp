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
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
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

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::SKALL]);
  xmlNodePtr node = doc.getRoot();
  UPtr<StructureFactorInput> sf_in;
  sf_in = std::make_unique<StructureFactorInput>(node);
  StructureFactorEstimator sfe(*sf_in, pset_ions, pset_elec);
}

TEST_CASE("StructureFactorEstimator::Accumulate", "[estimators]")
{
  using Input = qmcplusplus::testing::ValidStructureFactorInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::SKALL]);
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

  std::vector<ParticleSet::ParticlePos> deterministic_rs = {
      {
          {0.5156778693, 0.9486072659, -1.174232483},
          {-0.3166678548, 1.376550555, 1.445290089},
          {1.960713625, 2.472656727, 1.051449418},
          {0.7458532453, 0.5551358461, 4.080774784},
          {-0.3515016139, -0.5192222595, 0.9941511154},
          {-0.8354426622, 0.7071638107, -0.3409843445},
          {0.4386044741, 1.237378716, 2.331874132},
          {2.125850677, 0.322106719, 0.5825730562},
      },
      {
          {-0.4633736908, 0.06318772584, -0.8580153584},
          {-1.1749264, -0.6276503801, 0.07458759099},
          {1.327618122, 2.085829258, 1.415749788},
          {0.9114726782, 0.1789183617, -0.08135545254},
          {-2.267908812, 0.8029287457, 0.9522812963},
          {1.50271523, -1.844935298, 0.2805620432},
          {3.168934584, 0.1348338127, 1.371092796},
          {0.8310229182, 1.070827127, 1.180167317},
      },
      {
          {-0.04847732186, -1.201739907, -1.700527787},
          {0.1589259505, -0.3096047044, -2.06662631},
          {2.2559762, 1.62913239, -0.8024448156},
          {2.5347929, 3.121092796, 1.609780669},
          {-0.2892376184, -0.1520225108, -2.727613688},
          {0.2477154732, 0.5039232969, 2.995702744},
          {3.679345131, 3.037770271, 2.808899403},
          {0.6418578625, 1.935944557, 1.515637875},
      },
  };
  std::vector<ParticleSet> psets =
      testing::generateRandomParticleSets(pset_elec, pset_ions, deterministic_rs, nwalkers, generate_test_data);

  auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);

  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);

  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& spomap = wavefunction_pool.getWaveFunction("wavefunction")->getSPOMap();

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  // These are just empty arguments to hang the accumulation test, StructureFactorEstimator never accesses into them.
  // In the application the estimator manager calls accumulate and all these vectors are really populated.
  std::vector<QMCHamiltonian> hams;
  auto ref_wfns = convertUPtrToRefVector(twfcs);
  auto ref_hams = makeRefVector<QMCHamiltonian>(hams);

  FakeRandom<OHMMS_PRECISION_FULL> rng;

  auto updateWalker = [](auto& walker, auto& pset_target, auto& trial_wavefunction) {
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    //pset_target.saveWalker(walker);
  };

  for (int iw = 0; iw < nwalkers; ++iw)
    updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  // using QMCT = OperatorEstBase::QMCT;
  // std::vector<QMCT::RealType> rng_reals(nwalkers * QMCT::DIM * 2);

  // if (generate_test_data)
  // {
  //   testing::RandomForTest<QMCT::RealType> rft;
  //   rft.fillVecRng(rng_reals);
  //   std::cout << "rng_reals = " << NativePrint(rng_reals) << '\n';
  // }
  // else
  // {
  //   rng_reals = {
  //       0.6121701598,  0.120757781,  0.1690697521, 0.3017894626, 0.4360590279, 0.8539064527,
  //       0.7692624927,  0.2977159917, 0.295325309,  0.3778624535, 0.149162963,  0.2090236992,
  //       0.02247832529, 0.6495640278, 0.4202244878, 0.6017025113, 0.2386821359, 0.7122730613,
  //   };
  // }
  // auto it_rng_reals = rng_reals.begin();

  // for (ParticleSet& pset : psets)
  // {
  //   pset.R[0] = ParticleSet::PosType(*it_rng_reals++, *it_rng_reals++, *it_rng_reals++);
  //   pset.R[1] = ParticleSet::PosType(*it_rng_reals++, *it_rng_reals++, *it_rng_reals++);
  // }
  // for (int iw = 0; iw < nwalkers; ++iw)
  //   updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  CHECK(sfe.getNumKPoints() == 608);

  decltype(KContainer::kpts) kpoint_lists = {
      {-1, -1, -1}, {-1, 0, 0},   {0, -1, 0},   {0, 0, -1},   {0, 0, 1},    {0, 1, 0},    {1, 0, 0},    {1, 1, 1},
      {-1, -1, 0},  {-1, 0, -1},  {0, -1, -1},  {0, 1, 1},    {1, 0, 1},    {1, 1, 0},    {-2, -1, -1}, {-1, -2, -1},
      {-1, -1, -2}, {-1, 0, 1},   {-1, 1, 0},   {0, -1, 1},   {0, 1, -1},   {1, -1, 0},   {1, 0, -1},   {1, 1, 2},
      {1, 2, 1},    {2, 1, 1},    {-2, -2, -1}, {-2, -1, -2}, {-2, -1, 0},  {-2, 0, -1},  {-1, -2, -2}, {-1, -2, 0},
      {-1, -1, 1},  {-1, 0, -2},  {-1, 1, -1},  {-1, 1, 1},   {0, -2, -1},  {0, -1, -2},  {0, 1, 2},    {0, 2, 1},
      {1, -1, -1},  {1, -1, 1},   {1, 0, 2},    {1, 1, -1},   {1, 2, 0},    {1, 2, 2},    {2, 0, 1},    {2, 1, 0},
      {2, 1, 2},    {2, 2, 1},    {-2, -2, -2}, {-2, 0, 0},   {0, -2, 0},   {0, 0, -2},   {0, 0, 2},    {0, 2, 0},
      {2, 0, 0},    {2, 2, 2},    {-2, -2, 0},  {-2, 0, -2},  {0, -2, -2},  {0, 2, 2},    {2, 0, 2},    {2, 2, 0},
      {-3, -2, -2}, {-2, -3, -2}, {2, 3, 2},    {3, 2, 2},    {-3, -1, -1}, {-2, -2, -3}, {-2, 0, 1},   {-2, 1, 0},
      {-1, -3, -1}, {-1, -1, -3}, {-1, 0, 2},   {-1, 2, 0},   {0, -2, 1},   {0, -1, 2},   {0, 1, -2},   {0, 2, -1},
      {1, -2, 0},   {1, 0, -2},   {1, 1, 3},    {1, 3, 1},    {2, -1, 0},   {2, 0, -1},   {2, 2, 3},    {3, 1, 1},
      {-3, -2, -1}, {-3, -1, -2}, {-2, -3, -1}, {-2, -1, -3}, {-2, -1, 1},  {-2, 1, -1},  {-1, -3, -2}, {-1, -2, -3},
      {-1, -2, 1},  {-1, 1, -2},  {-1, 1, 2},   {-1, 2, 1},   {1, -2, -1},  {1, -1, -2},  {1, -1, 2},   {1, 2, -1},
      {1, 2, 3},    {1, 3, 2},    {2, -1, 1},   {2, 1, -1},   {2, 1, 3},    {2, 3, 1},    {3, 1, 2},    {3, 2, 1},
      {-2, -3, -3}, {2, 3, 3},    {-3, -3, -2}, {-3, -2, -3}, {-3, -1, 0},  {-3, 0, -1},  {-2, 1, 1},   {-1, -3, 0},
      {-1, -1, 2},  {-1, 0, -3},  {-1, 2, -1},  {0, -3, -1},  {0, -1, -3},  {0, 1, 3},    {0, 3, 1},    {1, -2, 1},
      {1, 0, 3},    {1, 1, -2},   {1, 3, 0},    {2, -1, -1},  {3, 0, 1},    {3, 1, 0},    {3, 2, 3},    {3, 3, 2},
      {-3, -3, -1}, {-3, -2, 0},  {-3, -1, -3}, {-3, 0, -2},  {-2, -3, 0},  {-2, -2, 1},  {-2, 0, -3},  {-2, 1, -2},
      {-1, -3, -3}, {-1, 2, 2},   {0, -3, -2},  {0, -2, -3},  {0, 2, 3},    {0, 3, 2},    {1, -2, -2},  {1, 3, 3},
      {2, -1, 2},   {2, 0, 3},    {2, 2, -1},   {2, 3, 0},    {3, 0, 2},    {3, 1, 3},    {3, 2, 0},    {3, 3, 1},
      {-3, -3, -3}, {-3, 0, 0},   {0, -3, 0},   {0, 0, -3},   {0, 0, 3},    {0, 3, 0},    {3, 0, 0},    {3, 3, 3},
      {-4, -2, -2}, {-2, -4, -2}, {-2, -2, -4}, {-2, 0, 2},   {-2, 2, 0},   {0, -2, 2},   {0, 2, -2},   {2, -2, 0},
      {2, 0, -2},   {2, 2, 4},    {2, 4, 2},    {4, 2, 2},    {-4, -2, -3}, {-4, -2, -1}, {-4, -1, -2}, {-3, -2, -4},
      {-3, -1, 1},  {-3, 1, -1},  {-2, -4, -3}, {-2, -4, -1}, {-2, -3, -4}, {-2, -1, -4}, {-2, -1, 2},  {-2, 1, 2},
      {-2, 2, -1},  {-2, 2, 1},   {-1, -4, -2}, {-1, -3, 1},  {-1, -2, -4}, {-1, -2, 2},  {-1, 1, -3},  {-1, 1, 3},
      {-1, 2, -2},  {-1, 3, 1},   {1, -3, -1},  {1, -2, 2},   {1, -1, -3},  {1, -1, 3},   {1, 2, -2},   {1, 2, 4},
      {1, 3, -1},   {1, 4, 2},    {2, -2, -1},  {2, -2, 1},   {2, -1, -2},  {2, 1, -2},   {2, 1, 4},    {2, 3, 4},
      {2, 4, 1},    {2, 4, 3},    {3, -1, 1},   {3, 1, -1},   {3, 2, 4},    {4, 1, 2},    {4, 2, 1},    {4, 2, 3},
      {-4, -3, -2}, {-3, -4, -2}, {3, 4, 2},    {4, 3, 2},    {-3, -4, -3}, {-3, -3, -4}, {-3, -3, 0},  {-3, 0, -3},
      {-3, 0, 1},   {-3, 1, 0},   {-1, -4, -1}, {-1, -1, -4}, {-1, 0, 3},   {-1, 3, 0},   {0, -3, -3},  {0, 3, 3},
      {1, -3, 0},   {1, 0, -3},   {1, 1, 4},    {1, 4, 1},    {3, -1, 0},   {3, 0, -1},   {3, 0, 3},    {3, 3, 0},
      {3, 3, 4},    {3, 4, 3},    {-4, -3, -3}, {-4, -1, -1}, {0, -3, 1},   {0, -1, 3},   {0, 1, -3},   {0, 3, -1},
      {4, 1, 1},    {4, 3, 3},    {-3, -2, 1},  {-2, -3, 1},  {2, 3, -1},   {3, 2, -1},   {-4, -1, -3}, {-3, -1, -4},
      {-3, 1, -2},  {-2, 1, -3},  {-1, -4, -3}, {-1, -3, -4}, {-1, 2, 3},   {-1, 3, 2},   {1, -3, -2},  {1, -2, -3},
      {1, 3, 4},    {1, 4, 3},    {2, -1, 3},   {3, -1, 2},   {3, 1, 4},    {4, 1, 3},    {-4, -3, -1}, {-3, -4, -1},
      {3, 4, 1},    {4, 3, 1},    {-4, -4, -3}, {-4, -3, -4}, {-4, -1, 0},  {-4, 0, -1},  {-3, -4, -4}, {-3, 1, 1},
      {-1, -4, 0},  {-1, -1, 3},  {-1, 0, -4},  {-1, 3, -1},  {0, -4, -1},  {0, -1, -4},  {0, 1, 4},    {0, 4, 1},
      {1, -3, 1},   {1, 0, 4},    {1, 1, -3},   {1, 4, 0},    {3, -1, -1},  {3, 4, 4},    {4, 0, 1},    {4, 1, 0},
      {4, 3, 4},    {4, 4, 3},    {-4, -4, -2}, {-4, -2, -4}, {-4, -2, 0},  {-4, 0, -2},  {-2, -4, -4}, {-2, -4, 0},
      {-2, -2, 2},  {-2, 0, -4},  {-2, 2, -2},  {-2, 2, 2},   {0, -4, -2},  {0, -2, -4},  {0, 2, 4},    {0, 4, 2},
      {2, -2, -2},  {2, -2, 2},   {2, 0, 4},    {2, 2, -2},   {2, 4, 0},    {2, 4, 4},    {4, 0, 2},    {4, 2, 0},
      {4, 2, 4},    {4, 4, 2},    {-4, -4, -4}, {-4, 0, 0},   {0, -4, 0},   {0, 0, -4},   {0, 0, 4},    {0, 4, 0},
      {4, 0, 0},    {4, 4, 4},    {-5, -3, -3}, {-3, -5, -3}, {-3, -3, -5}, {-3, 0, 2},   {-3, 2, 0},   {-2, -2, -5},
      {-2, 0, 3},   {-2, 3, 0},   {0, -3, 2},   {0, -2, 3},   {0, 2, -3},   {0, 3, -2},   {2, -3, 0},   {2, 0, -3},
      {2, 2, 5},    {3, -2, 0},   {3, 0, -2},   {3, 3, 5},    {3, 5, 3},    {5, 3, 3},    {-5, -2, -2}, {-2, -5, -2},
      {2, 5, 2},    {5, 2, 2},    {-4, -4, -1}, {-4, -3, 0},  {-4, -1, -4}, {-4, 0, -3},  {-3, -4, 0},  {-3, -3, 1},
      {-3, 0, -4},  {-3, 1, -3},  {-1, -4, -4}, {-1, 3, 3},   {0, -4, -3},  {0, -3, -4},  {0, 3, 4},    {0, 4, 3},
      {1, -3, -3},  {1, 4, 4},    {3, -1, 3},   {3, 0, 4},    {3, 3, -1},   {3, 4, 0},    {4, 0, 3},    {4, 1, 4},
      {4, 3, 0},    {4, 4, 1},    {-3, 2, -1},  {-2, -5, -3}, {-2, 3, 1},   {2, -3, -1},  {2, 5, 3},    {3, -2, 1},
      {-5, -3, -2}, {-5, -2, -3}, {-3, -5, -2}, {-3, -2, -5}, {-3, -1, 2},  {-2, -3, -5}, {-2, 1, 3},   {-1, -3, 2},
      {-1, 2, -3},  {1, -2, 3},   {1, 3, -2},   {2, -1, -3},  {2, 3, 5},    {3, 1, -2},   {3, 2, 5},    {3, 5, 2},
      {5, 2, 3},    {5, 3, 2},    {-4, -5, -3}, {-4, -1, 1},  {-4, 1, -1},  {-2, 3, -1},  {2, -3, 1},   {4, -1, 1},
      {4, 1, -1},   {4, 5, 3},    {-5, -4, -3}, {-3, -5, -4}, {-3, 2, 1},   {-1, -5, -2}, {-1, -4, 1},  {-1, 4, 1},
      {1, -4, -1},  {1, 4, -1},   {1, 5, 2},    {3, -2, -1},  {3, 5, 4},    {5, 4, 3},    {-5, -3, -4}, {-5, -1, -2},
      {-4, -3, -5}, {-3, -4, -5}, {-3, 1, 2},   {-2, -1, -5}, {-2, -1, 3},  {-1, -2, -5}, {-1, -2, 3},  {-1, 1, -4},
      {-1, 1, 4},   {-1, 3, -2},  {1, -3, 2},   {1, -1, -4},  {1, -1, 4},   {1, 2, -3},   {1, 2, 5},    {2, 1, -3},
      {2, 1, 5},    {3, -1, -2},  {3, 4, 5},    {4, 3, 5},    {5, 1, 2},    {5, 3, 4},    {-2, -5, -1}, {2, 5, 1},
      {-5, -2, -1}, {5, 2, 1},    {-5, -4, -4}, {-5, -1, -1}, {-4, -5, -4}, {-4, -4, -5}, {-4, 0, 1},   {-4, 1, 0},
      {-1, -5, -1}, {-1, -1, -5}, {-1, 0, 4},   {-1, 4, 0},   {0, -4, 1},   {0, -1, 4},   {0, 1, -4},   {0, 4, -1},
      {1, -4, 0},   {1, 0, -4},   {1, 1, 5},    {1, 5, 1},    {4, -1, 0},   {4, 0, -1},   {4, 4, 5},    {4, 5, 4},
      {5, 1, 1},    {5, 4, 4},    {-5, -4, -2}, {-5, -3, -1}, {-5, -1, -3}, {-4, -5, -2}, {-4, 1, -2},  {-3, -5, -1},
      {-3, 2, -2},  {-2, -4, 1},  {-2, -3, 2},  {-2, 1, -4},  {-2, 2, -3},  {-2, 3, 2},   {-1, -5, -3}, {-1, 4, 2},
      {1, -4, -2},  {1, 5, 3},    {2, -3, -2},  {2, -2, 3},   {2, -1, 4},   {2, 3, -2},   {2, 4, -1},   {3, -2, 2},
      {3, 5, 1},    {4, -1, 2},   {4, 5, 2},    {5, 1, 3},    {5, 3, 1},    {5, 4, 2},    {-5, -2, -4}, {-4, -2, -5},
      {-4, -2, 1},  {-3, -2, 2},  {-3, -1, -5}, {-2, -5, -4}, {-2, -4, -5}, {-2, 2, 3},   {-1, -3, -5}, {-1, 2, 4},
      {1, -2, -4},  {1, 3, 5},    {2, -2, -3},  {2, 4, 5},    {2, 5, 4},    {3, 1, 5},    {3, 2, -2},   {4, 2, -1},
      {4, 2, 5},    {5, 2, 4},    {-4, -4, 0},  {-4, 0, -4},  {0, -4, -4},  {0, 4, 4},    {4, 0, 4},    {4, 4, 0},
      {-5, -5, -3}, {5, 5, 3},    {-5, -3, -5}, {-5, -2, 0},  {-5, 0, -2},  {-3, 2, 2},   {-2, -5, 0},  {-2, -2, 3},
      {-2, 0, -5},  {-2, 3, -2},  {0, -5, -2},  {0, 5, 2},    {2, -3, 2},   {2, 0, 5},    {2, 2, -3},   {2, 5, 0},
      {3, -2, -2},  {5, 0, 2},    {5, 2, 0},    {5, 3, 5},    {0, -2, -5},  {0, 2, 5},    {-3, -5, -5}, {3, 5, 5},
      {-5, -5, -4}, {-5, -4, -5}, {-5, -1, 0},  {-5, 0, -1},  {-4, -5, -5}, {-4, 1, 1},   {-1, -5, -4}, {-1, -5, 0},
      {-1, -1, 4},  {-1, 0, -5},  {-1, 4, -1},  {0, -5, -1},  {0, -1, -5},  {0, 1, 5},    {0, 5, 1},    {1, -4, 1},
      {1, 0, 5},    {1, 1, -4},   {1, 5, 0},    {1, 5, 4},    {4, -1, -1},  {4, 5, 5},    {5, 0, 1},    {5, 1, 0},
      {5, 4, 5},    {5, 5, 4},    {-5, -4, -1}, {-5, -1, -4}, {-4, -5, -1}, {-4, -3, 1},  {-4, -1, -5}, {-4, 1, -3},
      {-3, -4, 1},  {-3, 1, -4},  {-1, -4, -5}, {-1, 3, 4},   {-1, 4, 3},   {1, -4, -3},  {1, -3, -4},  {1, 4, 5},
      {3, -1, 4},   {3, 4, -1},   {4, -1, 3},   {4, 1, 5},    {4, 3, -1},   {4, 5, 1},    {5, 1, 4},    {5, 4, 1},
  };
  //std::cout << "kpoint_lists = " << NativePrint(sfe.getKLists().kpts) << '\n';

  double tolerance = 0.1;
  {
    INFO("Checking kpoint_lists");
    auto check = checkVector(sfe.getKLists().kpts, kpoint_lists, true);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }
  
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
    FAIL_CHECK("Test always fails when generating new test reference data.");
  }
  auto sfk_e_e_expected = StructureFactorTests::getSKElecElec<Value>();
  auto rhok_e_expected  = StructureFactorTests::getRhoKElec<Value>();

  tolerance = 1e-6;
  if constexpr (std::is_same_v<RealAlias<QMCTraits::ValueType>, float>)
    tolerance = 2e-6;
  {
    INFO("In sfk_e_e test");
    auto check = checkVector(sfk_e_e_expected, sfk_e_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }
  {
    INFO("In rhok_e test");
    auto check = checkVector(rhok_e_expected, rhok_e, true, tolerance);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }
}


} // namespace qmcplusplus
