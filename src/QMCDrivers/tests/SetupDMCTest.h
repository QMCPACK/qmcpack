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

#include "QMCApp/tests/MinimalParticlePool.h"
#include "QMCApp/tests/MinimalWaveFunctionPool.h"
#include "QMCApp/tests/MinimalHamiltonianPool.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/tests/SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
class SetupDMCTest : public SetupPools
{
public:
  SetupDMCTest(int nranks = 4)
      : population(comm->size(),
                   particle_pool->getParticleSet("e"),
                   wavefunction_pool->getPrimary(),
                   hamiltonian_pool->getPrimary(),
                   comm->rank()),
        num_ranks(nranks),
        qmcdrv_input(3)
  {
    if (Concurrency::maxThreads<>() < 8)
      num_crowds = Concurrency::maxThreads<>();
  }

  DMCBatched operator()()
  {
    Libxml2Document doc;
    doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
    node = doc.getRoot();

    qmcdrv_input.readXML(node);
    dmcdrv_input.readXML(node);

    QMCDriverInput qmc_input_copy(qmcdrv_input);
    DMCDriverInput dmc_input_copy(dmcdrv_input);
    return {std::move(qmc_input_copy),
            std::move(dmc_input_copy),
            population,
            *(wavefunction_pool->getPrimary()),
            *(hamiltonian_pool->getPrimary()),
            *(wavefunction_pool),
            comm};
  }

public:
  MCPopulation population;

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
