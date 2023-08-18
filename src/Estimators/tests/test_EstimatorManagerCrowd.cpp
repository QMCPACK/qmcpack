//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/EstimatorManagerCrowd.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/tests/EstimatorManagerNewTest.h"
#include "Estimators/EstimatorInputDelegates.h"
#include "EstimatorManagerInputTest.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/ProjectData.h"

namespace qmcplusplus
{
TEST_CASE("EstimatorManagerCrowd::EstimatorManagerCrowd", "[estimators]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;

  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());


  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& pset            = *(particle_pool.getParticleSet("e"));
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  auto& twf             = *(wavefunction_pool.getWaveFunction("wavefunction"));
  auto& ham             = *(hamiltonian_pool.getPrimary());

  EstimatorManagerNew emn(comm, std::move(emi), ham, pset, twf);

  CHECK(emn.getNumEstimators() == 2);
  CHECK(emn.getNumScalarEstimators() == 0);

  EstimatorManagerCrowd emc(emn);
}

TEST_CASE("EstimatorManagerCrowd PerParticleHamiltonianLogger integration", "[estimators]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;

  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewVMCInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());


  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto& pset = *(particle_pool.getParticleSet("e"));
  // This is where the pset properties "properies" gain the different hamiltonian operator values.
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  auto& twf = *(wavefunction_pool.getWaveFunction("wavefunction"));
  auto& ham = *(hamiltonian_pool.getPrimary());

  ham.informOperatorsOfListener();

  PerParticleHamiltonianLoggerInput pphli;
  emi.append(std::move(pphli));

  EstimatorManagerNew emn(comm, std::move(emi), ham, pset, twf);

  CHECK(emn.getNumEstimators() == 3);
  CHECK(emn.getNumScalarEstimators() == 0);

  // We repeat a bunch of test_QMCHamiltonian here to get things set up.

  UPtrVector<QMCHamiltonian> hams;
  UPtrVector<TrialWaveFunction> twfs;
  std::vector<ParticleSet> psets;

  int num_walkers   = 4;
  int num_electrons = particle_pool.getParticleSet("e")->getTotalNum();
  int num_ions      = particle_pool.getParticleSet("ion")->getTotalNum();

  for (int iw = 0; iw < num_walkers; ++iw)
  {
    psets.emplace_back(pset);
    psets.back().randomizeFromSource(*particle_pool.getParticleSet("ion"));
    twfs.emplace_back(twf.makeClone(psets.back()));
    hams.emplace_back(hamiltonian_pool.getPrimary()->makeClone(psets.back(), *twfs.back()));
  }

  EstimatorManagerCrowd emc(emn);

  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  std::vector<MCPWalker> walkers(num_walkers, MCPWalker(pset.getTotalNum()));

  for (auto& walker : walkers)
  {
    walker.R          = pset.R;
    walker.spins      = pset.spins;
    walker.Properties = pset.Properties;
    walker.registerData();
    walker.DataSet.allocate();
  }

  auto walker_refs = makeRefVector<MCPWalker>(walkers);
  RefVectorWithLeader<MCPWalker> walker_list{walker_refs[0], walker_refs};

  auto p_refs = makeRefVector<ParticleSet>(psets);
  RefVectorWithLeader<ParticleSet> p_list{p_refs[0], p_refs};

  ResourceCollection pset_res("test_pset_res");
  p_list.getLeader().createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  auto twf_refs = convertUPtrToRefVector(twfs);
  RefVectorWithLeader<TrialWaveFunction> twf_list{twf_refs[0], twf_refs};

  ResourceCollection wfc_res("test_wfc_res");
  twf_list.getLeader().createResource(wfc_res);
  ResourceCollectionTeamLock<TrialWaveFunction> mw_wfc_lock(wfc_res, twf_list);

  auto ham_refs = convertUPtrToRefVector(hams);
  RefVectorWithLeader<QMCHamiltonian> ham_list(ham_refs[0], ham_refs);

  ResourceCollection ham_res("test_ham_res");
  ham_list.getLeader().createResource(ham_res);
  ResourceCollectionTeamLock<QMCHamiltonian> mw_ham_lock(ham_res, ham_list);

  emc.registerListeners(ham_list);

  //   Setup RNG
  FakeRandom<OHMMS_PRECISION_FULL> rng;

  // Without this QMCHamiltonian::mw_evaluate segfaults
  // Because the CoulombPBCAA hamiltonian component has PtclRhoK (StructFact) that is invalid.
  ParticleSet::mw_update(p_list);

  QMCHamiltonian::mw_evaluate(ham_list, twf_list, p_list);

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
    ham.saveProperty(walker.getPropertyBase());
  };
  for (int iw = 0; iw < num_walkers; ++iw)
    savePropertiesIntoWalker(*(hams[iw]), walkers[iw]);

  emc.accumulate(walker_refs, p_refs, twf_refs, rng);
}


} // namespace qmcplusplus
