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
  QMCDriverPools(const RuntimeOptions& runtime_options, Communicate* comm)
      : particle(MinimalParticlePool::make_diamondC_1x1x1(comm)),
        wavefunction(MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle)),
        hamiltonian(MinimalHamiltonianPool::make_hamWithEE(comm, particle, wavefunction))
  {}
  ParticleSetPool particle;
  WaveFunctionPool wavefunction;
  HamiltonianPool hamiltonian;
};

class QMCDriverInputTestWrapper : public QMCDriver
{
public:
  QMCDriverInputTestWrapper(const ProjectData& project_data, QMCDriverPools& pools, Communicate* comm)
      : QMCDriver(project_data, *pools.particle.getWalkerSet("e"), pools.wavefunction.getWaveFunction().value(),
                  pools.hamiltonian.getHamiltonian().value(), comm, "QMCDriverInputTestWrapper")
  {}

  bool readXML(xmlNodePtr node) { return putQMCInfo(node); }
  IndexType getBlocksBetweenRecompute() const { return nBlocksBetweenRecompute; }

  void run() override {}
  bool put(xmlNodePtr) override { return true; }
  QMCRunType getRunType() override { return QMCRunType::VMC; }
};

auto createDriver(const RuntimeOptions& runtime_options,
                  Communicate* comm,
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

TEST_CASE("QMCDriver retired check-properties parameters", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm = OHMMS::Controller;
  ProjectData test_project("testing", ProjectData::DriverVersion::LEGACY);
  QMCDriverPools pools(test_project.getRuntimeOptions(), comm);

#ifdef MIXED_PRECISION
  constexpr int default_blocks_between_recompute = 1;
#else
  constexpr int default_blocks_between_recompute = 10;
#endif

  SECTION("omitted")
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(R"(<qmc method="vmc" move="pbyp"/>)"));
    QMCDriverInputTestWrapper driver(test_project, pools, comm);
    REQUIRE(driver.readXML(doc.getRoot()));
    CHECK(driver.getBlocksBetweenRecompute() == default_blocks_between_recompute);
  }

  for (const std::string& alias : {"checkProperties", "checkproperties", "check_properties"})
  {
    CAPTURE(alias);
    const std::string input =
        "<qmc method=\"vmc\" move=\"pbyp\">"
        "<parameter name=\"" +
        alias +
        "\">not-an-integer</parameter>"
        "<parameter name=\"blocks_between_recompute\">3</parameter>"
        "</qmc>";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input));
    QMCDriverInputTestWrapper driver(test_project, pools, comm);
    REQUIRE(driver.readXML(doc.getRoot()));
    CHECK(driver.getBlocksBetweenRecompute() == 3);
  }
}

TEST_CASE("QMCDriverFactory retired mode attributes", "[qmcapp]")
{
  ProjectData test_project("testing", ProjectData::DriverVersion::LEGACY);
  QMCDriverFactory driver_factory(test_project);

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(
      R"(<qmc method="vmc" move="pbyp" multiple="yes" warp="yes" append="yes" profiling="yes"/>)"));
  const auto das = driver_factory.readSection(doc.getRoot());

  CHECK(das.new_run_type == QMCRunType::VMC);
  CHECK(das.what_to_do[UPDATE_MODE]);
  CHECK(das.append_run);
  CHECK(das.enable_profiling);
}

TEST_CASE("QMCDriverFactory create VMC Driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  using DV = ProjectData::DriverVersion;
  ProjectData test_project("testing", DV::LEGACY);
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
  Communicate* comm;
  comm = OHMMS::Controller;
  using namespace testing;

  SECTION("driver version behavior")
  {
    ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
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

TEST_CASE("QMCDriverFactory create DMC driver", "[qmcapp]")
{
  using namespace testing;
  Communicate* comm;
  comm = OHMMS::Controller;

  using DV = ProjectData::DriverVersion;
  ProjectData test_project("testing", DV::LEGACY);
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
  Communicate* comm;
  comm = OHMMS::Controller;

  SECTION("driver version behavior")
  {
    ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
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
