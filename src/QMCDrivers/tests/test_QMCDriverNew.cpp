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


#include "catch.hpp"

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/tests/QMCDriverNewTestWrapper.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "Message/Communicate.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "EstimatorInputDelegates.h"
#include "QMCDrivers/MCPopulation.h"
#include "Concurrency/Info.hpp"
#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
TEST_CASE("QMCDriverNew tiny case", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_tiny_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmcdriver(test_project, std::move(qmcdriver_input), walker_confs,
                                    MCPopulation(comm->size(), comm->rank(), particle_pool.getParticleSet("e"),
                                                 wavefunction_pool.getPrimary(), hamiltonian_pool.getPrimary()),
                                    comm);

  // setStatus must be called before process
  std::string root_name{"Test"};
  //For later sections this appears to contain important state.
  std::string prev_config_file{""};

  qmcdriver.setStatus(root_name, prev_config_file, false);
  // We want to express out expectations of the QMCDriver state machine so we catch
  // changes to it over time.
  outputManager.resume();

  REQUIRE(qmcdriver.getBranchEngine() == nullptr);
  qmcdriver.process(node);
  REQUIRE(qmcdriver.get_num_living_walkers() == 1);

  // What else should we expect after process
}

#ifdef _OPENMP
TEST_CASE("QMCDriverNew more crowds than threads", "[drivers]")
{
  using namespace testing;

  Concurrency::OverrideMaxCapacity<> override(8);
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  int num_crowds = 9;

  // test is a no op except for openmp, max threads is >> than num cores
  // in other concurrency models.
  if (Concurrency::maxCapacity<>() != 8)
    throw std::runtime_error("Insufficient threads available to match test input");

  QMCDriverInput qmcdriver_copy(qmcdriver_input);
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmc_batched(test_project, std::move(qmcdriver_copy), walker_confs,
                                      MCPopulation(comm->size(), comm->rank(), particle_pool.getParticleSet("e"),
                                                   wavefunction_pool.getPrimary(), hamiltonian_pool.getPrimary()),
                                      comm);
  QMCDriverNewTestWrapper::TestNumCrowdsVsNumThreads<ParallelExecutor<>> testNumCrowds;
  testNumCrowds(9);
  testNumCrowds(8);
}

TEST_CASE("QMCDriverNew walker counts", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_dmc_input_sections[valid_dmc_input_dmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);

  int num_crowds = 8;

  if (Concurrency::maxCapacity<>() < 8)
    num_crowds = Concurrency::maxCapacity<>();

  if (num_crowds < 8)
    throw std::runtime_error("Insufficient threads available to match test input");

  QMCDriverInput qmcdriver_copy(qmcdriver_input);
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmc_batched(test_project, std::move(qmcdriver_copy), walker_confs,
                                      MCPopulation(comm->size(), comm->rank(), particle_pool.getParticleSet("e"),
                                                   wavefunction_pool.getPrimary(), hamiltonian_pool.getPrimary()),
                                      comm);

  qmc_batched.testAdjustGlobalWalkerCount();
}
#endif

TEST_CASE("QMCDriverNew test driver operations", "[drivers]")
{
  using namespace testing;
  Concurrency::OverrideMaxCapacity<> override(8);
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* comm = OHMMS::Controller;
  outputManager.pause();

  Libxml2Document doc;
  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_tiny_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(node);
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);

  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  WalkerConfigurations walker_confs;
  QMCDriverNewTestWrapper qmcdriver(test_project, std::move(qmcdriver_input), walker_confs,
                                    MCPopulation(comm->size(), comm->rank(), particle_pool.getParticleSet("e"),
                                                 wavefunction_pool.getPrimary(), hamiltonian_pool.getPrimary()),
                                    comm);


  auto tau       = 1.0;
  auto invmass   = 0.2;
  auto spin_mass = 0.5;
  {
    TauParams<QMCTraits::RealType, CoordsType::POS> taus(tau, invmass, spin_mass);
    constexpr auto mct = CoordsType::POS;
    auto mc_coords     = MCCoords<mct>(3);
    QMCTraits::PosType p1({0.1, 0.2, 0.3});
    QMCTraits::PosType p2({0.4, 0.5, 0.6});
    QMCTraits::PosType p3({0.7, 0.8, 0.9});
    mc_coords.positions = {p1, p2, p3};

    auto deltas = mc_coords;
    QMCDriverNew::scaleBySqrtTau(taus, deltas);
    mc_coords += deltas;
    CHECK(Approx(mc_coords.positions[0][0]) == 0.14472135955);
    CHECK(Approx(mc_coords.positions[0][1]) == 0.28944271910);
    CHECK(Approx(mc_coords.positions[0][2]) == 0.43416407865);
    CHECK(Approx(mc_coords.positions[1][0]) == 0.57888543820);
    CHECK(Approx(mc_coords.positions[1][1]) == 0.72360679775);
    CHECK(Approx(mc_coords.positions[1][2]) == 0.86832815730);
    CHECK(Approx(mc_coords.positions[2][0]) == 1.01304951685);
    CHECK(Approx(mc_coords.positions[2][1]) == 1.15777087640);
    CHECK(Approx(mc_coords.positions[2][2]) == 1.30249223595);

    std::vector<QMCTraits::RealType> loggf(3), loggb(3);

    qmcdriver.computeLogGreensFunction(deltas, taus, loggf);
    CHECK(Approx(loggf[0]) == -0.07);
    CHECK(Approx(loggf[1]) == -0.385);
    CHECK(Approx(loggf[2]) == -0.97);

    qmcdriver.computeLogGreensFunction(mc_coords, taus, loggb);
    CHECK(Approx(loggb[0]) == -0.733049516850);
    CHECK(Approx(loggb[1]) == -4.031772342675);
    CHECK(Approx(loggb[2]) == -10.15797187635);
  }

  {
    TauParams<QMCTraits::RealType, CoordsType::POS_SPIN> taus(tau, invmass, spin_mass);
    constexpr auto mct = CoordsType::POS_SPIN;
    auto mc_coords     = MCCoords<mct>(3);
    QMCTraits::PosType p1({-0.1, -0.2, -0.3});
    QMCTraits::PosType p2({-0.4, -0.5, -0.6});
    QMCTraits::PosType p3({-0.7, -0.8, -0.9});
    mc_coords.positions = {p1, p2, p3};
    mc_coords.spins     = {0.1, 0.2, 0.3};

    auto deltas = mc_coords;
    QMCDriverNew::scaleBySqrtTau(taus, deltas);
    mc_coords += deltas;
    CHECK(Approx(mc_coords.positions[0][0]) == -0.14472135955);
    CHECK(Approx(mc_coords.positions[0][1]) == -0.28944271910);
    CHECK(Approx(mc_coords.positions[0][2]) == -0.43416407865);
    CHECK(Approx(mc_coords.positions[1][0]) == -0.57888543820);
    CHECK(Approx(mc_coords.positions[1][1]) == -0.72360679775);
    CHECK(Approx(mc_coords.positions[1][2]) == -0.86832815730);
    CHECK(Approx(mc_coords.positions[2][0]) == -1.01304951685);
    CHECK(Approx(mc_coords.positions[2][1]) == -1.15777087640);
    CHECK(Approx(mc_coords.positions[2][2]) == -1.30249223595);

    CHECK(Approx(mc_coords.spins[0]) == 0.163245553203);
    CHECK(Approx(mc_coords.spins[1]) == 0.326491106407);
    CHECK(Approx(mc_coords.spins[2]) == 0.489736659610);

    std::vector<QMCTraits::RealType> loggf(3), loggb(3);

    qmcdriver.computeLogGreensFunction(deltas, taus, loggf);
    CHECK(Approx(loggf[0]) == -0.075);
    CHECK(Approx(loggf[1]) == -0.405);
    CHECK(Approx(loggf[2]) == -1.015);

    qmcdriver.computeLogGreensFunction(mc_coords, taus, loggb);
    CHECK(Approx(loggb[0]) == -0.766360905151);
    CHECK(Approx(loggb[1]) == -4.165017895878);
    CHECK(Approx(loggb[2]) == -10.457774371057);
  }

  {
    outputManager.resume();
    qmcdriver.testMeasureImbalance();
    outputManager.pause();
  }
}

} // namespace qmcplusplus
