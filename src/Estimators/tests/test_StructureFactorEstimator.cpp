//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "StructureFactorEstimator.h"
#include "ValidStructureFactorInput.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "GenerateRandomParticleSets.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include <iostream>

namespace qmcplusplus
{

constexpr bool generate_test_data = true;

namespace testing
{
class StructureFactorAccess
{
  using Real = StructureFactorEstimator::Real;
public:
  const Vector<Real>& getSKElecElec(const StructureFactorEstimator& sfe) const { return sfe.getSKElecElec(); }
  const Vector<std::complex<Real>>& getRhoKElec(const StructureFactorEstimator& sfe) const { return sfe.getRhoKElec(); }
};
} // namespace testing

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

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  StructureFactorEstimator sfe(*sf_in, pset_ions, pset_elec);

  std::vector<OperatorEstBase::MCPWalker> walkers;
  int nwalkers = 3;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);
    
  std::vector<ParticleSet::ParticlePos> deterministic_rs = {{
                                                                {-0.6759092808, 0.835668385, 1.985307097},
                                                                {0.09710352868, -0.76751858, -1.89306891},
                                                                {-0.5605484247, -0.9578875303, 1.476860642},
                                                                {2.585144997, 1.862680197, 3.282609463},
                                                                {-0.1961335093, 1.111888766, -0.578481257},
                                                                {1.794641614, 1.6000278, -0.9474347234},
                                                                {2.157717228, 0.9254754186, 2.263158321},
                                                                {1.883366346, 2.136350632, 3.188981533},
                                                            },
                                                            {
                                                                {-0.2079261839, -0.2796236873, 0.5512072444},
                                                                {-0.2823159397, 0.7537326217, 0.01526880637},
                                                                {3.533515453, 2.433290243, 0.9281452894},
                                                                {2.051767349, 2.312927485, 0.7089259624},
                                                                {-1.043096781, 0.8190526962, -0.1958218962},
                                                                {0.9210210443, 0.7726522088, 0.3962054551},
                                                                {2.043324947, 0.3482068777, 3.39059639},
                                                                {0.9103830457, 2.167978764, 2.341906071},
                                                            },
                                                            {
                                                                {-0.466550231, 0.09173964709, -0.3779250085},
                                                                {-0.4211375415, -2.017466068, -1.691870451},
                                                                {2.090800285, 1.88529861, 2.152359247},
                                                                {2.973145723, 1.718174577, 3.822324753},
                                                                {-0.8552014828, -0.3484517336, -0.2870049179},
                                                                {0.2349359095, -0.5025780797, 0.2305756211},
                                                                {-0.03547382355, 2.279159069, 3.057915211},
                                                                {2.535993099, 1.637133598, 3.689830303},
                                                            }};

  std::vector<ParticleSet> psets =
      testing::generateRandomParticleSets(pset_ions, pset_elec, deterministic_rs, nwalkers, generate_test_data);

  auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);

  // These are just empty arguments to hang the accumulation test, StructureFactorEstimator never accesses into them.
  // In the application the estimator manager calls accumulate and all these vectors are really populated.
  std::vector<TrialWaveFunction> wfns;
  std::vector<QMCHamiltonian> hams;
  auto ref_wfns    = makeRefVector<TrialWaveFunction>(wfns);
  auto ref_hams    = makeRefVector<QMCHamiltonian>(hams);

  FakeRandom<OHMMS_PRECISION_FULL> rng;

  sfe.accumulate(ref_walkers, ref_psets, ref_wfns, ref_hams, rng);

  testing::StructureFactorAccess sfa;  
  auto& sfk_e_e = sfa.getSKElecElec(sfe);
  auto& rhok_e = sfa.getRhoKElec(sfe);

  if constexpr (generate_test_data) {
    std::cout << "sfk_e_e_exected = ";
    std::cout << NativePrint(sfk_e_e) << '\n';
    std::cout << "rhok_e_expected = ";
    std::cout << NativePrint(rhok_e) << '\n';
    FAIL_CHECK("Test always fails when generating new test reference data.");
  }
}

} // namespace qmcplusplus
