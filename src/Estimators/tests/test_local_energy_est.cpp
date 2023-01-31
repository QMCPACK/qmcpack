//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "QMCDrivers/WalkerProperties.h"
#include "io/hdf/hdf_archive.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

TEST_CASE("LocalEnergyOnly", "[estimators]")
{
  LocalEnergyOnlyEstimator le_est;

  const SimulationCell simulation_cell;
  MCWalkerConfiguration W(simulation_cell);
  W.setName("electrons");
  W.create({1});
  W.createWalkers(1);

  (*W.begin())->Properties(WP::LOCALENERGY) = 1.1;

  le_est.accumulate(W, W.begin(), W.end(), 1.0);

  CHECK(le_est.getName() == "LocalEnergyOnlyEstimator");
  CHECK(le_est.scalars[0].mean() == Approx(1.1));
}

TEST_CASE("LocalEnergy", "[estimators]")
{
  QMCHamiltonian H;
  LocalEnergyEstimator le_est(H, false);

  CHECK(le_est.getName() == "LocalEnergyEstimator");

  std::unique_ptr<LocalEnergyEstimator> le_est2{le_est.clone()};
  REQUIRE(le_est2 != nullptr);
  REQUIRE(le_est2.get() != &le_est);

  const SimulationCell simulation_cell;
  MCWalkerConfiguration W(simulation_cell);
  W.setName("electrons");
  W.create({1});
  W.createWalkers(1);

  (*W.begin())->Properties(WP::LOCALENERGY)    = 1.1;
  (*W.begin())->Properties(WP::LOCALPOTENTIAL) = 1.2;

  le_est.accumulate(W, W.begin(), W.end(), 1.0);

  // from enum in LocalEnergyEstimator
  // 0 - ENERGY_INDEX
  // 1 - ENERGY2_INDEX
  // 2 - POTENTIAL_INDEX
  CHECK(le_est.scalars[0].mean() == Approx(1.1));
  REQUIRE(le_est.scalars[1].mean() == le_est.scalars[0].mean2());
  CHECK(le_est.scalars[2].mean() == Approx(1.2));
}

TEST_CASE("LocalEnergy with hdf5", "[estimators]")
{
  QMCHamiltonian H;
  LocalEnergyEstimator le_est(H, true);

  const SimulationCell simulation_cell;
  MCWalkerConfiguration W(simulation_cell);
  W.setName("electrons");
  W.create({1});
  W.createWalkers(1);

  (*W.begin())->Properties(WP::LOCALENERGY)    = 1.1;
  (*W.begin())->Properties(WP::LOCALPOTENTIAL) = 1.2;

  std::vector<ObservableHelper> h5desc;

  std::filesystem::path filename("tmp_obs.h5");
  hdf_archive h_file;
  h_file.create(filename);
  le_est.registerObservables(h5desc, h_file);
  h_file.close();
  REQUIRE(std::filesystem::exists(filename));
  // Check contents?
  REQUIRE(std::filesystem::remove(filename));

  LocalEnergyEstimator::RecordListType record;
  le_est.add2Record(record);
  // Not sure how to test this - for now make sure it doesn't crash
}

} // namespace qmcplusplus
