//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "QMCDrivers/QMCDriverFactory.h"
#include "DMC/DMCBatched.h"
#include "ValidQMCInputSections.h"
#include "tests/EstimatorManagerInputTest.h"
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
#include <EstimatorInputDelegates.h>
#include "DMCBatchedTestAccessor.h"
namespace qmcplusplus
{

TEST_CASE("DMCBatched::estimator_measurement_period", "[drivers]")
{
  using namespace testing;
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();
  using DMCInput = qmcplusplus::testing::DmcInputs;
  Libxml2Document doc;
  bool okay = doc.parseFromString(DMCInput::getXml(DMCInput::valid::EMPERIOD));
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  QMCDriverFactory driver_factory(test_project);
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  auto qmc_system = particle_pool.getWalkerSet("e");

  // This code repeats code in driver factory because it really only
  // exists because of legacy and erases driver type.
  QMCDriverInput qmcdriver_input;
  DMCDriverInput dmcdriver_input;
  qmcdriver_input.readXML(cur);
  dmcdriver_input.readXML(cur);

  DMCBatched(test_project,
             std::move(qmcdriver_input,
                       makeEstimatorManager(std::nullopt, qmcdriver_input.get_estimator_manager_input(), qmc_system,
                                            *primaryPsi, *primaryH, particle_pool.getPool(), comm), ),
             std::move(dmcdriver_input), *qmc_system,
             MCPopulation(comm->size(), comm->rank(), &qmc_system, primaryPsi, primaryH),
             RandomNumberControl::getChildrenRefs(), comm);

  qmc_driver->process(node);

  qmc_driver->setStatus(test_project.currentMainRoot(), "", false);

  using Dmcbta = testing::DMCBatchedTestAccessor;
  Dmcbta::mockRunStart(*qmc_driver);
}

} // namespace qmcplusplus
