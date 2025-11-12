//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "EnergyDensityTest.h"
#include "Utilities/ProjectData.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/StdRandom.h"
#include "GenerateRandomParticleSets.h"
#include "ValidEnergyDensityInput.h"

namespace qmcplusplus
{
namespace testing
{
std::function<MockGoldWalkerElements(Communicate*, RuntimeOptions)> make_gold_walker_elem_ee =
    makeGoldWalkerElementsWithEE;

EnergyDensityTest::EnergyDensityTest(Communicate* comm, int num_walkers, bool generate_test_data)
    : EnergyDensityTest(comm, num_walkers, make_gold_walker_elem_ee, generate_test_data)
{}

EnergyDensityTest::EnergyDensityTest(Communicate* comm,
                                     int num_walkers,
                                     std::function<MockGoldWalkerElements(Communicate*, RuntimeOptions)> make_gold_elem,
                                     bool generate_test_data)
    : test_project_("test", ProjectData::DriverVersion::BATCH),
      gold_elem_(make_gold_elem(comm, test_project_.getRuntimeOptions())),
      walkers_(num_walkers, MCPWalker(gold_elem_.pset_elec.getTotalNum()))
{
  // This has to be first because before the golden Hamiltonian is copied
  // it needs to know if there are per particle listeners.
  gold_elem_.ham.informOperatorsOfListener();

  psets_ = testing::generateRandomParticleSets(gold_elem_.pset_elec, gold_elem_.pset_ions, deterministic_rs_,
                                               num_walkers, generate_test_data);

  auto pset_refs(makeRefVector<ParticleSet>(psets_));
  auto& trial_wavefunction = gold_elem_.twf;

  Libxml2Document doc;
  using Input = testing::EnergyDensityInputs;
  Input input;
  bool okay       = doc.parseFromString(Input::getXml(Input::valid::CELL));
  xmlNodePtr node = doc.getRoot();

  UPtr<EnergyDensityInput> edein;
  edein = std::make_unique<EnergyDensityInput>(node);
  NEEnergyDensityEstimator e_den_est(*edein, gold_elem_.particle_pool.getPool());

  for (int iw = 0; iw < num_walkers; ++iw)
  {
    twfs_.emplace_back(gold_elem_.twf.makeClone(psets_[iw]));
    hams_.emplace_back(gold_elem_.ham.makeClone(psets_[iw], *twfs_.back()));
  }

  auto walker_refs = makeRefVector<MCPWalker>(walkers_);

  RefVector<QMCHamiltonian> ham_refs = convertUPtrToRefVector(hams_);

  hams_[0]->createResource(ham_res_);
  psets_[0].createResource(pset_res_);
  twfs_[0]->createResource(wfc_res_);
}

NEEnergyDensityEstimator& EnergyDensityTest::getEnergyDensityEstimator()
{
  if (!eden_est_)
  {
    Libxml2Document doc;
    using Input     = testing::EnergyDensityInputs;
    bool okay       = doc.parseFromString(Input::getXml(Input::valid::CELL));
    xmlNodePtr node = doc.getRoot();
    edein_          = std::make_unique<EnergyDensityInput>(node);
    eden_est_       = std::make_unique<NEEnergyDensityEstimator>(*edein_, gold_elem_.particle_pool.getPool());
  }
  return *eden_est_;
}

RefVectorWithLeader<EnergyDensityTest::MCPWalker> EnergyDensityTest::getWalkerList()
{
  return {walkers_[0], makeRefVector<MCPWalker>(walkers_)};
}
RefVectorWithLeader<ParticleSet> EnergyDensityTest::getPSetList()
{
  return {psets_[0], makeRefVector<ParticleSet>(psets_)};
}
RefVectorWithLeader<QMCHamiltonian> EnergyDensityTest::getHamList()
{
  return {*hams_[0], convertUPtrToRefVector(hams_)};
}
RefVectorWithLeader<TrialWaveFunction> EnergyDensityTest::getTwfList()
{
  return {*twfs_[0], convertUPtrToRefVector(twfs_)};
}

ResourceCollection& EnergyDensityTest::getPSetRes() { return pset_res_; }
ResourceCollection& EnergyDensityTest::getHamRes() { return ham_res_; }
ResourceCollection& EnergyDensityTest::getTwfRes() { return wfc_res_; }


} // namespace testing
} // namespace qmcplusplus
