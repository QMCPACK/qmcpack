//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OneBodyDensityMatrices.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "TrialWaveFunction.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/StdRandom.h"
#include "Utilities/StlPrettyPrint.hpp"
#include <iostream>

namespace qmcplusplus
{

namespace testing
{
template<typename T>
class OneBodyDensityMatricesTests
{
public:
  using Evaluators  = OneBodyDensityMatricesInput::Evaluator;
  using Integrators = OneBodyDensityMatricesInput::Integrator;
  using Sampling    = OneBodyDensityMatrices::Sampling;
  using MCPWalker   = OneBodyDensityMatrices::MCPWalker;

  OneBodyDensityMatricesTests() = default;
  void testGenerateSamples(onebodydensitymatrices::Inputs input,
                           OneBodyDensityMatrices& obdm,
                           ParticleSet& pset_target,
                           StdRandom<T>& rng)
  {
    using namespace onebodydensitymatrices;
    switch (input)
    {
    case (valid_obdm_input):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 64);
      break;
    case (valid_obdm_input_scale):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 0);
      break;
    case (valid_obdm_input_grid):
      obdm.generateSamples(1.0, pset_target, rng);
      CHECK(obdm.nmoves_ == 0);
      CHECK(obdm.samples_ == pow(22, OHMMS_DIM));
      break;
    }
  }

  void testEvaluateMatrix(OneBodyDensityMatrices& obdm,
                          ParticleSet& pset,
                          TrialWaveFunction& trial_wavefunction,
                          MCPWalker& walker,
                          StdRandom<T>& rng)
  {
    obdm.evaluateMatrix(pset, trial_wavefunction, walker, rng);
  }

  void dumpData(OneBodyDensityMatrices& obdm)
  {
    std::cout << *(obdm.data_);
  }
};

} // namespace testing


TEST_CASE("OneBodyDensityMatrices::OneBodyDensityMatrices", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  auto lattice     = testing::makeTestLattice();
  auto species_set = testing::makeSpeciesSet(SpeciesCases::GOOD);

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  {
    // Good constructor
    OneBodyDensityMatrices obDenMat(std::move(obdmi), lattice, species_set, wf_factory, pset_target);
    // Good copy constructor
    OneBodyDensityMatrices obDenMat2(obDenMat);
  }
  {
    species_set = testing::makeSpeciesSet(SpeciesCases::NO_MEMBERSIZE);
    CHECK_THROWS_AS(OneBodyDensityMatrices(std::move(obdmi), lattice, species_set, wf_factory, pset_target),
                    UniformCommunicateError);
  }

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::generateSamples", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;

  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  auto& wf_factory  = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  auto samplingCaseRunner = [&pset_target, &species_set, &wf_factory](Inputs test_case) {
    Libxml2Document doc;

  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[test_case]);
  if (!okay)
      throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);

    OneBodyDensityMatrices obDenMat(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);

    OneBodyDensityMatricesTests<double> obdmt;
    //Get control over which rng is used.
    //we don't want FakeRandom.
    StdRandom<double> rng;
    obdmt.testGenerateSamples(test_case, obDenMat, pset_target, rng);
  };

  samplingCaseRunner(valid_obdm_input);
  samplingCaseRunner(valid_obdm_input_scale);
  samplingCaseRunner(valid_obdm_input_grid);

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::clone()", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;

  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  auto& wf_factory  = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[Inputs::valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);

  OneBodyDensityMatrices original(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);
  auto clone = original.clone();
  REQUIRE(clone != nullptr);
  REQUIRE(clone.get() != &original);
  REQUIRE(dynamic_cast<decltype(&original)>(clone.get()) != nullptr);
}

TEST_CASE("OneBodyDensityMatrices::accumulate", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  using MCPWalker = OperatorEstBase::MCPWalker;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  OneBodyDensityMatrices(std::move(obdmi), pset_target.Lattice, species_set, wf_factory, pset_target);

  std::vector<MCPWalker> walkers;
  int nwalkers = 4;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  std::vector<ParticleSet> psets;
  for (int iw = 0; iw < nwalkers; ++iw)
    psets.emplace_back(pset_target);

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> twfcs(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    twfcs[iw] = trial_wavefunction.makeClone(psets[iw]);

  // now the framework for testing accumulation is done

  outputManager.resume();
}

TEST_CASE("OneBodyDensityMatrices::evaluateMatrix", "[estimators]")
{
  using namespace testing;
  using namespace onebodydensitymatrices;
  using MCPWalker = OperatorEstBase::MCPWalker;
  using namespace onebodydensitymatrices;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[valid_obdm_input]);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset        = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset.getSpeciesSet();
  OneBodyDensityMatrices obdm(std::move(obdmi), pset.Lattice, species_set, wf_factory);
  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  StdRandom<double> rng;
  MCPWalker walker;
  walker.Weight = 0.95;
  OneBodyDensityMatricesTests<double> obdmt;
  trial_wavefunction.evaluateLog(pset);
  obdmt.testEvaluateMatrix(obdm, pset, trial_wavefunction, walker, rng);
  obdmt.dumpData(obdm);
  // now the framework for testing accumulation is done
  outputManager.resume();
}

} // namespace qmcplusplus
