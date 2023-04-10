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
#include "EstimatorInputDelegates.h"
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
  VMCBatchedTest(const ProjectData& project_data) : project_data_(project_data)
  {
    Concurrency::OverrideMaxCapacity<> override(8);
    comm_ = OHMMS::Controller;
  }

  void testCalcDefaultLocalWalkers()
  {
    using namespace testing;
    Concurrency::OverrideMaxCapacity<> override(8);

    Libxml2Document doc;
    doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]);
    xmlNodePtr node = doc.getRoot();
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(node);

    auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm_);
    auto wavefunction_pool =
        MinimalWaveFunctionPool::make_diamondC_1x1x1(project_data_.getRuntimeOptions(), comm_, particle_pool);
    auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm_, particle_pool, wavefunction_pool);
  }

private:
  Communicate* comm_;
  const ProjectData& project_data_;
};
} // namespace testing

TEST_CASE("VMCBatched::calc_default_local_walkers", "[drivers]")
{
  using namespace testing;
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  VMCBatchedTest vbt(test_project);
  vbt.testCalcDefaultLocalWalkers();
}

} // namespace qmcplusplus
