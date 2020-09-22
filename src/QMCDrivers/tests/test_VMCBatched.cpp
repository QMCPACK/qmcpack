//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Message/Communicate.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"
#include "Particle/SampleStack.h"

namespace qmcplusplus
{
namespace testing
{
class VMCBatchedTest
{
public:
  VMCBatchedTest()
  {
    Concurrency::OverrideMaxCapacity<> override(8);
    Communicate* comm;
    OHMMS::Controller->initialize(0, NULL);
    comm_ = OHMMS::Controller;
  }

  void testCalcDefaultLocalWalkers()
  {
    using namespace testing;
    Concurrency::OverrideMaxCapacity<> override(8);
    Communicate* comm;
    OHMMS::Controller->initialize(0, NULL);
    comm = OHMMS::Controller;

    Libxml2Document doc;
    doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]);
    xmlNodePtr node = doc.getRoot();
    QMCDriverInput qmcdriver_input(3);
    qmcdriver_input.readXML(node);

    MinimalParticlePool mpp;
    ParticleSetPool particle_pool = mpp(comm);
    MinimalWaveFunctionPool wfp;
    WaveFunctionPool wavefunction_pool = wfp(comm, &particle_pool);
    wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
    MinimalHamiltonianPool mhp;
    HamiltonianPool hamiltonian_pool = mhp(comm, &particle_pool, &wavefunction_pool);


      
  }

private:
  Communicate* comm_;
};
} // namespace testing

TEST_CASE("VMCBatched::calc_default_local_walkers", "[drivers]")
{
  using namespace testing;
  VMCBatchedTest vbt;
  vbt.testCalcDefaultLocalWalkers();
}

} // namespace qmcplusplus
