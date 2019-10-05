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

#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCApp/QMCDriverFactory.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include "Message/Communicate.h"

namespace qmcplusplus
{
TEST_CASE("QMCDriverFactory create VMC Driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  QMCDriverFactory driver_factory;

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCRunType::VMC);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm,&particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm,&particle_pool, &wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  std::unique_ptr<QMCDriverInterface> last_driver;
  std::unique_ptr<QMCDriverInterface> qmc_driver;

  MCPopulation population(comm->size(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());

  qmc_driver = driver_factory.newQMCDriver(std::move(last_driver), 0, node, das, *qmc_system, particle_pool,
                                           wavefunction_pool, hamiltonian_pool, population, comm);
  REQUIRE(qmc_driver != nullptr);
}

TEST_CASE("QMCDriverFactory create VMCBatched driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  QMCDriverFactory driver_factory;

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCRunType::VMC_BATCH);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  std::unique_ptr<QMCDriverInterface> last_driver;
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  MCPopulation population(comm->size(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());

  qmc_driver = driver_factory.newQMCDriver(std::move(last_driver), 0, node, das, *qmc_system, particle_pool,
                                           wavefunction_pool, hamiltonian_pool, population, comm);
  REQUIRE(qmc_driver != nullptr);
}

TEST_CASE("QMCDriverFactory create DMC driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  QMCDriverFactory driver_factory;

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCRunType::DMC);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  std::unique_ptr<QMCDriverInterface> last_driver;
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  MCPopulation population(comm->size(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());

  qmc_driver = driver_factory.newQMCDriver(std::move(last_driver), 0, node, das, *qmc_system, particle_pool,
                                           wavefunction_pool, hamiltonian_pool, population, comm);
  REQUIRE(qmc_driver != nullptr);
}

TEST_CASE("QMCDriverFactory create DMCBatched driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  QMCDriverFactory driver_factory;

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCRunType::DMC_BATCH);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  std::unique_ptr<QMCDriverInterface> last_driver;
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  MCPopulation population(comm->size(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());

  qmc_driver = driver_factory.newQMCDriver(std::move(last_driver), 0, node, das, *qmc_system, particle_pool,
                                           wavefunction_pool, hamiltonian_pool, population, comm);
  REQUIRE(qmc_driver != nullptr);
}

} // namespace qmcplusplus
