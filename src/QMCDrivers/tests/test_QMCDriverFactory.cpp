//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "ProjectData.h"
#include "Message/Communicate.h"
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/DMC/DMC.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/DMC/DMCBatched.h"

namespace qmcplusplus
{
namespace testing
{
class QMCDriverPools
{
public:
  QMCDriverPools(Communicate* comm)
      : particle(MinimalParticlePool::make_diamondC_1x1x1(comm)),
        wavefunction(MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle)),
        hamiltonian(MinimalHamiltonianPool::make_hamWithEE(comm, particle, wavefunction))
  {}
  ParticleSetPool particle;
  WaveFunctionPool wavefunction;
  HamiltonianPool hamiltonian;
};

auto createDriver(Communicate* comm,
                  QMCDriverFactory& driver_factory,
                  xmlNodePtr node,
                  QMCDriverFactory::DriverAssemblyState& das)
{
  QMCDriverPools dr_pools(comm);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = dr_pools.particle.getWalkerSet(target);
  return driver_factory.createQMCDriver(node, das, *qmc_system, dr_pools.particle, dr_pools.wavefunction,
                                        dr_pools.hamiltonian, comm);
}

} // namespace testing

TEST_CASE("QMCDriverFactory create VMC Driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  ProjectData test_project;
  QMCDriverFactory driver_factory(test_project);

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  REQUIRE(das.new_run_type == QMCRunType::VMC);

  // testing::QMCDriverPools dr_pools(comm);
  // std::string target("e");
  // MCWalkerConfiguration* qmc_system = dr_pools.particle.getWalkerSet(target);

  // std::unique_ptr<QMCDriverInterface> qmc_driver;

  // qmc_driver =
  //   driver_factory.createQMCDriver(node, das, *qmc_system, dr_pools.particle, dr_pools.wavefunction, dr_pools.hamiltonian, comm);

  auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);

  REQUIRE(qmc_driver != nullptr);
  REQUIRE_THROWS(dynamic_cast<VMCBatched&>(*qmc_driver));
  REQUIRE_NOTHROW(dynamic_cast<VMC&>(*qmc_driver));
  CHECK(qmc_driver->getEngineName() == "VMC");
}

TEST_CASE("QMCDriverFactory create VMCBatched driver", "[qmcapp]")
{
  Communicate* comm;
  comm = OHMMS::Controller;
  using namespace testing;

  SECTION("driver version behavior")
  {
    ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]);
    REQUIRE(okay);
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::VMC_BATCH);

    auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);
    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<VMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<VMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "VMCBatched");
  }
  SECTION("Deprecated _batch behavior")
  {
    using namespace testing;
    ProjectData test_project;
    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_batch_input_vmc_batch_index]);
    REQUIRE(okay);
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::VMC_BATCH);

    auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);

    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<VMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<VMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "VMCBatched");
  }
}

TEST_CASE("QMCDriverFactory create DMC driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  ProjectData test_project;
  QMCDriverFactory driver_factory(test_project);

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_index]);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  REQUIRE(das.new_run_type == QMCRunType::DMC);

  auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);

  REQUIRE(qmc_driver != nullptr);
  REQUIRE_THROWS(dynamic_cast<DMCBatched&>(*qmc_driver));
  REQUIRE_NOTHROW(dynamic_cast<DMC&>(*qmc_driver));
  CHECK(qmc_driver->getEngineName() == "DMC");
}

TEST_CASE("QMCDriverFactory create DMCBatched driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  SECTION("driver version behavior")
  {
    ProjectData test_project("test", ProjectData::DriverVersion::BATCH);

    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
    REQUIRE(okay);
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::DMC_BATCH);

    auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);

    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<DMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<DMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "DMCBatched");
  }
  SECTION("Deprecated _batch behavior")
  {
    ProjectData test_project;
    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_batch_input_dmc_batch_index]);
    REQUIRE(okay);
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::DMC_BATCH);

    auto qmc_driver = testing::createDriver(comm, driver_factory, node, das);

    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<DMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<DMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "DMCBatched");
  }
}

} // namespace qmcplusplus
