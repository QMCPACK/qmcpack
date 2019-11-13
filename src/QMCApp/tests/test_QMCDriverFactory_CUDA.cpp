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

#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCApp/QMCDriverFactory.h"
#include "QMCDrivers/QMCDriverInterface.h"
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
TEST_CASE("QMCDriverFactory create VMC_CUDA Driver", "[qmcapp]")
{
  Communicate* comm;
  OHMMS::Controller->initialize(0, NULL);
  comm = OHMMS::Controller;

  QMCDriverFactory driver_factory;
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
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCRunType::VMC);

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
  MinimalHamiltonianPool mhp;
  HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);
  std::string target("e");
  MCWalkerConfiguration* qmc_system = particle_pool.getWalkerSet(target);

  MCPopulation population(comm->size(), particle_pool.getParticleSet("e"), wavefunction_pool.getPrimary(),
                            hamiltonian_pool.getPrimary());

  std::unique_ptr<QMCDriverInterface> last_driver;
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  qmc_driver = driver_factory.newQMCDriver(std::move(last_driver), 0, node, das, *qmc_system, particle_pool,
                                           wavefunction_pool, hamiltonian_pool, population, comm);
  REQUIRE(qmc_driver != nullptr);
}

} // namespace qmcplusplus
