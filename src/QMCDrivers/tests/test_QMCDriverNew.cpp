//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/tests/QMCDriverNewTestWrapper.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "Message/Communicate.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "QMCDrivers/MCPopulation.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
TEST_CASE("QMCDriverNew tiny case", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_tiny_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  SampleStack samples;
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmcdriver(std::move(qmcdriver_input),
                                    MCPopulation(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"),
                                                 wavefunction_pool.getPrimary(),
                                                 wavefunction_pool.getWaveFunctionFactory("wavefunction"),
                                                 hamiltonian_pool.getPrimary()),
                                    samples, comm);

  // setStatus must be called before process
  std::string root_name{"Test"};
  //For later sections this appears to contain important state.
  std::string prev_config_file{""};

  qmcdriver.setStatus(root_name, prev_config_file, false);
  // We want to express out expectations of the QMCDriver state machine so we catch
  // changes to it over time.
  outputManager.resume();

  REQUIRE(qmcdriver.getBranchEngine() == nullptr);
  qmcdriver.process(node);
  REQUIRE(qmcdriver.get_num_living_walkers() == 1);

  // What else should we expect after process
}

#ifdef _OPENMP
TEST_CASE("QMCDriverNew more crowds than threads", "[drivers]")
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
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  int num_crowds = 9;

  // test is a no op except for openmp, max threads is >> than num cores
  // in other concurrency models.
  if (Concurrency::maxCapacity<>() != 8)
    throw std::runtime_error("Insufficient threads available to match test input");

  QMCDriverInput qmcdriver_copy(qmcdriver_input);
  SampleStack samples;
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmc_batched(std::move(qmcdriver_copy),
                                      MCPopulation(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"),
                                                   wavefunction_pool.getPrimary(),
                                                   wavefunction_pool.getWaveFunctionFactory("wavefunction"),
                                                   hamiltonian_pool.getPrimary()),
                                      samples, comm);
  QMCDriverNewTestWrapper::TestNumCrowdsVsNumThreads<ParallelExecutor<>> testNumCrowds;
  testNumCrowds(9);
  testNumCrowds(8);
}

TEST_CASE("QMCDriverNew walker counts", "[drivers]")
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
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  int num_crowds = 8;

  if (Concurrency::maxCapacity<>() < 8)
    num_crowds = Concurrency::maxCapacity<>();

  if (num_crowds < 8)
    throw std::runtime_error("Insufficient threads available to match test input");

  QMCDriverInput qmcdriver_copy(qmcdriver_input);
  SampleStack samples;
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmc_batched(std::move(qmcdriver_copy),
                                      MCPopulation(1, comm->rank(), walker_confs, particle_pool.getParticleSet("e"),
                                                   wavefunction_pool.getPrimary(),
                                                   wavefunction_pool.getWaveFunctionFactory("wavefunction"),
                                                   hamiltonian_pool.getPrimary()),
                                      samples, comm);

  qmc_batched.testAdjustGlobalWalkerCount();
}
#endif

} // namespace qmcplusplus
