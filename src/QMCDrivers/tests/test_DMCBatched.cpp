//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch.hpp>

#include "Message/Communicate.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCDrivers/tests/SetupDMCTest.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "Platforms/Host/OutputManager.h"

namespace qmcplusplus
{
namespace testing
{
class DMCBatchedTest
{
public:
  DMCBatchedTest() { up_dtest_ = std::make_unique<SetupDMCTest>(1); }

  void testDependentObjectsValidAfterPopulationChange()
  {
    using namespace testing;
    SetupDMCTest& dtest = get_dtest();
  }

  SetupDMCTest& get_dtest() { return *up_dtest_; }

private:
  UPtr<SetupDMCTest> up_dtest_;
};
} // namespace testing

/** Since we check the DMC only feature of reserve walkers perhaps this should be
 *  a DMC integration test.
 */
TEST_CASE("DMCDriver+QMCDriverNew integration", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input(3);
  qmcdriver_input.readXML(node);
  DMCDriverInput dmcdriver_input;
  dmcdriver_input.readXML(node);
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));

  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, particle_pool, wavefunction_pool);
  SampleStack samples;
  WalkerConfigurations walker_confs;
  DMCBatched dmcdriver(std::move(qmcdriver_input), std::move(dmcdriver_input),
                       MCPopulation(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"),
                                    wavefunction_pool.getPrimary(), hamiltonian_pool.getPrimary())


                           ,
                       *(wavefunction_pool.getPrimary()), *(hamiltonian_pool.getPrimary()), comm);

  // setStatus must be called before process
  std::string root_name{"Test"};
  //For later sections this appears to contain important state.
  std::string prev_config_file{""};

  dmcdriver.setStatus(root_name, prev_config_file, false);
  // We want to express out expectations of the QMCDriver state machine so we catch
  // changes to it over time.
  CHECK(dmcdriver.getNewBranchEngine() == nullptr);
  outputManager.resume();

  dmcdriver.process(node);
  CHECK(dmcdriver.getNewBranchEngine() != nullptr);
  CHECK(dmcdriver.get_num_living_walkers() == 8);
  const QMCTraits::IndexType reserved_walkers = dmcdriver.get_num_living_walkers() + dmcdriver.get_num_dead_walkers();
  CHECK(reserved_walkers == 10);
  // What else should we expect after process
}


} // namespace qmcplusplus
