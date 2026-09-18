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
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/Crowd.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "EstimatorInputDelegates.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/StdRandom.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include "Utilities/for_testing/RandomForTest.h"
#include <MinimalHamiltonianPool.h>
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include "Particle/SampleStack.h"

#include <cmath>
#include <limits>

namespace qmcplusplus
{
namespace testing
{
class VMCBatchedTest
{
public:
  static void runAllParticleNoDrift()
  {
    constexpr bool generate_test_data = false;
    constexpr const char* allp_input  = R"(
      <qmc method="vmc_batch" move="allp" checkpoint="-1">
        <parameter name="crowds">1</parameter>
        <parameter name="total_walkers">2</parameter>
        <parameter name="warmupsteps">0</parameter>
        <parameter name="substeps">2</parameter>
        <parameter name="steps">1</parameter>
        <parameter name="blocks">1</parameter>
        <parameter name="timestep">0.01</parameter>
        <parameter name="usedrift">no</parameter>
      </qmc>
    )";

    ProjectData project_data("test", ProjectData::DriverVersion::BATCH);
    Communicate* comm = OHMMS::Controller;
    RandomForTest<QMCTraits::FullPrecRealType> random_for_test;
    Random.init(static_cast<int>(random_for_test() * 1000000));

    Libxml2Document doc;
    REQUIRE(doc.parseFromString(allp_input));
    xmlNodePtr node = doc.getRoot();

    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(node);
    VMCDriverInput vmcdriver_input;
    vmcdriver_input.readXML(node);
    REQUIRE(qmcdriver_input.get_update_mode() == "allp");
    REQUIRE_FALSE(qmcdriver_input.areWalkersSerialized());
    REQUIRE_FALSE(vmcdriver_input.get_use_drift());

    auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
    auto wavefunction_pool =
        MinimalWaveFunctionPool::make_diamondC_1x1x1(project_data.getRuntimeOptions(), comm, particle_pool);
    auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

    using RNG = RandomBase<QMCTraits::FullPrecRealType>;
    StdRandom<QMCTraits::FullPrecRealType> rng0(
        static_cast<StdRandom<QMCTraits::FullPrecRealType>::uint_type>(random_for_test() * 1000000));
    RefVector<RNG> rng_refs;
    rng_refs.push_back(rng0);

    WalkerConfigurations walker_confs;
    SampleStack samples;
    VMCBatched vmc_driver(project_data, std::move(qmcdriver_input), nullptr, std::move(vmcdriver_input), walker_confs,
                          MCPopulation(comm->size(), comm->rank(), *particle_pool.getParticleSet("e"),
                                       wavefunction_pool.getWaveFunction().value(),
                                       hamiltonian_pool.getHamiltonian().value()),
                          rng_refs, samples, comm);

    vmc_driver.setStatus("vmcbatched_allp", "", false);
    vmc_driver.process(node);

    REQUIRE(vmc_driver.get_num_living_walkers() == 2);
    REQUIRE(vmc_driver.crowds_.size() == 1);
    CHECK(vmc_driver.crowds_[0]->size() == 2);

    vmc_driver.run();
    CHECK(walker_confs.getActiveWalkers() == 2);
    std::vector<QMCTraits::FullPrecRealType> final_log_psis;
    std::vector<int> move_accepts;
    std::vector<int> move_rejects;
    for (const auto& crowd : vmc_driver.crowds_)
    {
      const int particles_per_walker = crowd->get_walker_elecs()[0].get().getTotalNum();
      CHECK(crowd->get_accept() + crowd->get_reject() == 2 * crowd->size() * particles_per_walker);
      move_accepts.push_back(crowd->get_accept());
      move_rejects.push_back(crowd->get_reject());
      for (const auto& twf : crowd->get_walker_twfs())
      {
        CHECK(std::isfinite(twf.get().getLogPsi()));
        final_log_psis.push_back(twf.get().getLogPsi());
      }
    }

    if constexpr (generate_test_data)
      app_log() << "VMCBatched allp move_accepts = " << NativePrint<std::vector<int>>(move_accepts)
                << "\nVMCBatched allp move_rejects = " << NativePrint<std::vector<int>>(move_rejects)
                << "\nVMCBatched allp final_log_psis = "
                << NativePrint<std::vector<QMCTraits::FullPrecRealType>>(final_log_psis) << '\n';
    else
    {
      const std::vector<int> expected_move_accepts{16};
      const std::vector<int> expected_move_rejects{16};
      const std::vector<QMCTraits::FullPrecRealType> expected_final_log_psis{-4.2715065656, -3.6563663443};
      const auto log_psi_epsilon = std::numeric_limits<QMCTraits::FullPrecRealType>::epsilon() * 1000000;

      CHECK(move_accepts == expected_move_accepts);
      CHECK(move_rejects == expected_move_rejects);
      REQUIRE(final_log_psis.size() == expected_final_log_psis.size());
      for (size_t iw = 0; iw < final_log_psis.size(); ++iw)
        CHECK(final_log_psis[iw] == Approx(expected_final_log_psis[iw]).epsilon(log_psi_epsilon));
    }
  }

  static void testAllParticleModeValidation()
  {
    constexpr const char* unsupported_allp_input = R"(
      <qmc method="vmc_batch" move="allp">
        <parameter name="usedrift">yes</parameter>
      </qmc>
    )";

    Libxml2Document doc;
    REQUIRE(doc.parseFromString(unsupported_allp_input));
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(doc.getRoot());
    VMCDriverInput vmcdriver_input;
    vmcdriver_input.readXML(doc.getRoot());

    REQUIRE_THROWS_WITH(VMCBatched::validateAllParticleMode(qmcdriver_input, vmcdriver_input, false),
                        "VMCBatched all-particle moves do not support drift.");

    constexpr const char* serialized_allp_input = R"(
      <qmc method="vmc_batch" move="allp">
        <parameter name="crowd_serialize_walkers">yes</parameter>
        <parameter name="usedrift">no</parameter>
      </qmc>
    )";
    Libxml2Document serialized_doc;
    REQUIRE(serialized_doc.parseFromString(serialized_allp_input));
    qmcdriver_input.readXML(serialized_doc.getRoot());
    vmcdriver_input.readXML(serialized_doc.getRoot());
    REQUIRE_THROWS_WITH(VMCBatched::validateAllParticleMode(qmcdriver_input, vmcdriver_input, false),
                        "VMCBatched all-particle moves do not support serialized crowd walkers.");

    constexpr const char* position_allp_input = R"(
      <qmc method="vmc_batch" move="allp">
        <parameter name="usedrift">no</parameter>
      </qmc>
    )";
    Libxml2Document position_doc;
    REQUIRE(position_doc.parseFromString(position_allp_input));
    QMCDriverInput position_qmcdriver_input;
    position_qmcdriver_input.readXML(position_doc.getRoot());
    VMCDriverInput position_vmcdriver_input;
    position_vmcdriver_input.readXML(position_doc.getRoot());
    REQUIRE_THROWS_WITH(VMCBatched::validateAllParticleMode(position_qmcdriver_input, position_vmcdriver_input, true),
                        "VMCBatched all-particle moves support position coordinates only.");

    constexpr const char* pbyp_input = R"(
      <qmc method="vmc_batch" move="pbyp">
        <parameter name="crowd_serialize_walkers">yes</parameter>
        <parameter name="usedrift">yes</parameter>
      </qmc>
    )";
    Libxml2Document pbyp_doc;
    REQUIRE(pbyp_doc.parseFromString(pbyp_input));
    QMCDriverInput pbyp_qmcdriver_input;
    pbyp_qmcdriver_input.readXML(pbyp_doc.getRoot());
    VMCDriverInput pbyp_vmcdriver_input;
    pbyp_vmcdriver_input.readXML(pbyp_doc.getRoot());
    REQUIRE_NOTHROW(VMCBatched::validateAllParticleMode(pbyp_qmcdriver_input, pbyp_vmcdriver_input, true));
  }
};
} // namespace testing

TEST_CASE("VMCBatched all-particle moves without drift", "[drivers]")
{
  testing::VMCBatchedTest::runAllParticleNoDrift();
}

TEST_CASE("VMCBatched all-particle mode validation", "[drivers]")
{
  testing::VMCBatchedTest::testAllParticleModeValidation();
}

} // namespace qmcplusplus
