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

#ifndef QMCPLUSPLUS_SETUP_DMCTEST_H
#define QMCPLUSPLUS_SETUP_DMCTEST_H

#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/SetupPools.h"
#include "ProjectData.h"
#include "tests/ValidQMCInputSections.h"

namespace qmcplusplus
{
namespace testing
{
class SetupDMCTest : public SetupPools
{
public:
  SetupDMCTest(int nranks = 4) : rng_pool(Concurrency::maxCapacity<>()), num_ranks(nranks), qmcdrv_input()
  {
    if (Concurrency::maxCapacity<>() < 8)
      num_crowds = Concurrency::maxCapacity<>();
  }

  DMCBatched operator()()
  {
    Libxml2Document doc;
    doc.parseFromString(DmcInputs::getXml(DmcInputs::valid::CROWDS));
    node = doc.getRoot();

    qmcdrv_input.readXML(node);
    dmcdrv_input.readXML(node);

    QMCDriverInput qmc_input_copy(qmcdrv_input);
    DMCDriverInput dmc_input_copy(dmcdrv_input);
    return {test_project,
            std::move(qmc_input_copy),
            nullptr,
            std::move(dmc_input_copy),
            walker_confs,
            MCPopulation(comm->size(), comm->rank(), particle_pool->getParticleSet("e"),
                         wavefunction_pool->getPrimary(), hamiltonian_pool->getPrimary()),
            rng_pool.getRngRefs(),
            comm};
  }

private:
  ProjectData test_project;

  RandomNumberGeneratorPool rng_pool;

public:
  WalkerConfigurations walker_confs;

  Libxml2Document doc;
  xmlNodePtr node;

  int num_ranks;
  int num_crowds = 8;

  QMCDriverInput qmcdrv_input;
  DMCDriverInput dmcdrv_input;
};
} // namespace testing
} // namespace qmcplusplus

#endif
