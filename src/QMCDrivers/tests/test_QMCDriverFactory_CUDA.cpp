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

/** @file
 *  Temporary DriverRactory test for legacy CUDA.
 * 
 *  Since there is not any plan to have the new batch driver and the legacy CUDA variant co-compile.
 *  This avoids polluting the new test code and new batch drivers with \#ifdef QMC_CUDA
 *  The new driver sources are simply not part of the build at all for the legacy CUDA
 *  variant.
 */

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
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "ProjectData.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
TEST_CASE("QMCDriverFactory create VMC_CUDA Driver", "[qmcapp]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  ProjectData test_project;
  QMCDriverFactory driver_factory(test_project);
  // clang-format off
  const char* driver_xml = R"(
  <qmc method="vmc" move="pbyp" gpu="yes">
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                1 </parameter>
    <parameter name="stepsbetweensamples">    1 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)";
  // clang-format on

  Libxml2Document doc;
  bool okay = doc.parseFromString(driver_xml);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(node);
  REQUIRE(das.new_run_type == QMCRunType::VMC);

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto hamiltonian_pool  = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  std::unique_ptr<QMCDriverInterface> qmc_driver;
  qmc_driver =
      driver_factory.createQMCDriver(node, das, *qmc_system, particle_pool, wavefunction_pool, hamiltonian_pool, comm);
  REQUIRE(qmc_driver != nullptr);
}

} // namespace qmcplusplus
