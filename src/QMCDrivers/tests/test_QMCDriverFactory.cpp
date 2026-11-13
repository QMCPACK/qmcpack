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
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/QMCDriverFactory.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
#include "ProjectData.h"
#include "Message/Communicate.h"
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/DMC/DMC.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "EstimatorInputDelegates.h"

namespace qmcplusplus
{
namespace testing
{
class QMCDriverPools
{
public:
  QMCDriverPools(const RuntimeOptions& runtime_options, Communicate& comm)
      : particle(MinimalParticlePool::make_diamondC_1x1x1(comm)),
        wavefunction(MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle)),
        hamiltonian(MinimalHamiltonianPool::make_hamWithEE(comm, particle, wavefunction))
  {}
  ParticleSetPool particle;
  WaveFunctionPool wavefunction;
  HamiltonianPool hamiltonian;
};

auto createDriver(const RuntimeOptions& runtime_options,
                  Communicate& comm,
                  QMCDriverFactory& driver_factory,
                  xmlNodePtr node,
                  QMCDriverFactory::DriverAssemblyState& das)
{
  QMCDriverPools dr_pools(runtime_options, comm);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = dr_pools.particle.getWalkerSet(target);
  return driver_factory.createQMCDriver(node, das, std::nullopt, *qmc_system, dr_pools.particle, dr_pools.wavefunction,
                                        dr_pools.hamiltonian, comm);
}

} // namespace testing

TEST_CASE("QMCDriverFactory create VMC Driver", "[qmcapp]")
{
  using namespace testing;
  Communicate& comm(*OHMMS::Controller);

  using DV = ProjectData::DriverVersion;
  ProjectData test_project(DV::LEGACY);
  QMCDriverFactory driver_factory(test_project);

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_index]));
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  REQUIRE(das.new_run_type == QMCRunType::VMC);

  auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);

  REQUIRE(qmc_driver != nullptr);
  REQUIRE_THROWS(dynamic_cast<VMCBatched&>(*qmc_driver));
  REQUIRE_NOTHROW(dynamic_cast<VMC&>(*qmc_driver));
  CHECK(qmc_driver->getEngineName() == "VMC");
}

TEST_CASE("QMCDriverFactory create VMCBatched driver", "[qmcapp]")
{
  Communicate& comm(*OHMMS::Controller);
  using namespace testing;

  SECTION("driver version behavior")
  {
    ProjectData test_project;
    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    REQUIRE(doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]));
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::VMC_BATCH);

    auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);
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
    REQUIRE(doc.parseFromString(valid_vmc_input_sections[valid_vmc_batch_input_vmc_batch_index]));
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::VMC_BATCH);

    auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);

    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<VMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<VMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "VMCBatched");
  }
}

TEST_CASE("QMCDriverFactory all-particle method support", "[qmcapp]")
{
  ProjectData project;
  QMCDriverFactory factory(project);
  Libxml2Document doc;

  SECTION("batched DMC rejects all-particle moves")
  {
    REQUIRE(doc.parseFromString("<qmc method=\"dmc\" move=\"allp\"/>"));
    REQUIRE_THROWS(factory.readSection(doc.getRoot()));
  }
  SECTION("batched optimization rejects all-particle moves")
  {
    REQUIRE(doc.parseFromString("<qmc method=\"linear\" move=\"allp\"/>"));
    REQUIRE_THROWS(factory.readSection(doc.getRoot()));
  }
  SECTION("legacy DMC retains all-particle support")
  {
    ProjectData legacy_project(ProjectData::DriverVersion::LEGACY);
    QMCDriverFactory legacy_factory(legacy_project);
    REQUIRE(doc.parseFromString("<qmc method=\"dmc\" move=\"allp\"/>"));
    REQUIRE_NOTHROW(legacy_factory.readSection(doc.getRoot()));
  }
}

TEST_CASE("QMCDriverFactory validates batched VMC all-particle options before construction", "[qmcapp]")
{
  using namespace testing;
  Communicate& comm(*OHMMS::Controller);
  ProjectData project;
  QMCDriverFactory factory(project);

  QMCDriverPools pools(project.getRuntimeOptions(), comm);
  MCWalkerConfiguration* qmc_system = pools.particle.getWalkerSet("e");
  const auto create_from_xml = [&](const char* xml, bool spinor = false) {
    qmc_system->setSpinor(spinor);
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(xml));
    xmlNodePtr node = doc.getRoot();
    auto das        = factory.readSection(node);
    return factory.createQMCDriver(node, das, std::nullopt, *qmc_system, pools.particle, pools.wavefunction,
                                   pools.hamiltonian, comm);
  };
  const auto require_error = [&](const char* xml, const char* expected_message, bool spinor = false) {
    try
    {
      create_from_xml(xml, spinor);
      FAIL("Expected factory validation to throw");
    }
    catch (const UniformCommunicateError& error)
    {
      CHECK(std::string(error.what()) == expected_message);
    }
  };

  SECTION("supported all-particle mode constructs batched VMC")
  {
    auto driver = create_from_xml("<qmc method=\"vmc_batch\" move=\"allp\"><parameter name=\"usedrift\">no</parameter></qmc>");
    REQUIRE(driver != nullptr);
    CHECK(driver->getEngineName() == "VMCBatched");
  }
  SECTION("drift is rejected")
  {
    require_error("<qmc method=\"vmc_batch\" move=\"allp\"><parameter name=\"usedrift\">yes</parameter></qmc>",
                  "VMCBatched all-particle moves do not support drift.");
  }
  SECTION("serialized walkers are rejected through the parsed public input")
  {
    require_error("<qmc method=\"vmc_batch\" move=\"allp\"><parameter name=\"usedrift\">no</parameter><parameter "
                  "name=\"crowd_serialize_walkers\">yes</parameter></qmc>",
                  "VMCBatched all-particle moves do not support serialized crowd walkers.");
  }
  SECTION("spinor particle systems are rejected")
  {
    require_error("<qmc method=\"vmc_batch\" move=\"allp\"><parameter name=\"usedrift\">no</parameter></qmc>",
                  "VMCBatched all-particle moves support position coordinates only.", true);
  }
}

TEST_CASE("QMCDriverFactory create DMC driver", "[qmcapp]")
{
  using namespace testing;
  Communicate& comm(*OHMMS::Controller);

  using DV = ProjectData::DriverVersion;
  ProjectData test_project(DV::LEGACY);
  QMCDriverFactory driver_factory(test_project);

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_index]));
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  REQUIRE(das.new_run_type == QMCRunType::DMC);

  auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);

  REQUIRE(qmc_driver != nullptr);
  REQUIRE_THROWS(dynamic_cast<DMCBatched&>(*qmc_driver));
  REQUIRE_NOTHROW(dynamic_cast<DMC&>(*qmc_driver));
  CHECK(qmc_driver->getEngineName() == "DMC");
}

TEST_CASE("QMCDriverFactory create DMCBatched driver", "[qmcapp]")
{
  using namespace testing;
  Communicate& comm(*OHMMS::Controller);

  SECTION("driver version behavior")
  {
    ProjectData test_project;
    QMCDriverFactory driver_factory(test_project);

    Libxml2Document doc;
    REQUIRE(doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]));
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::DMC_BATCH);

    auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);
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
    REQUIRE(doc.parseFromString(valid_dmc_input_sections[valid_dmc_batch_input_dmc_batch_index]));
    xmlNodePtr node                           = doc.getRoot();
    QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
    REQUIRE(das.new_run_type == QMCRunType::DMC_BATCH);

    auto qmc_driver = testing::createDriver(test_project.getRuntimeOptions(), comm, driver_factory, node, das);

    REQUIRE(qmc_driver != nullptr);
    REQUIRE_NOTHROW(dynamic_cast<DMCBatched&>(*qmc_driver));
    REQUIRE_THROWS(dynamic_cast<DMC&>(*qmc_driver));
    CHECK(qmc_driver->getEngineName() == "DMCBatched");
  }
}

} // namespace qmcplusplus
