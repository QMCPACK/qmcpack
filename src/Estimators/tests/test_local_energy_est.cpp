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


#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

TEST_CASE("LocalEnergyOnly", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  LocalEnergyOnlyEstimator le_est;

  MCWalkerConfiguration W;
  W.setName("electrons");
  W.create(1);
  W.createWalkers(1);

  (*W.begin())->Properties(LOCALENERGY) = 1.1;

  le_est.accumulate(W, W.begin(), W.end(), 1.0);

  REQUIRE(le_est.scalars[0].mean() == Approx(1.1));
}

TEST_CASE("LocalEnergy", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;


  QMCHamiltonian H;
  LocalEnergyEstimator le_est(H, false);

  LocalEnergyEstimator *le_est2 = dynamic_cast<LocalEnergyEstimator *>(le_est.clone());
  REQUIRE(le_est2 != NULL);
  REQUIRE(le_est2 != &le_est);

  MCWalkerConfiguration W;
  W.setName("electrons");
  W.create(1);
  W.createWalkers(1);

  (*W.begin())->Properties(LOCALENERGY) = 1.1;
  (*W.begin())->Properties(LOCALPOTENTIAL) = 1.2;

  le_est.accumulate(W, W.begin(), W.end(), 1.0);

  // from enum in LocalEnergyEstimator
  // 0 - ENERGY_INDEX
  // 1 - ENERGY2_INDEX
  // 2 - POTENTIAL_INDEX
  REQUIRE(le_est.scalars[0].mean() == Approx(1.1));
  REQUIRE(le_est.scalars[1].mean() == le_est.scalars[0].mean2());
  REQUIRE(le_est.scalars[2].mean() == Approx(1.2));
}

TEST_CASE("LocalEnergy with hdf5", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;


  QMCHamiltonian H;
  LocalEnergyEstimator le_est(H, true);

  MCWalkerConfiguration W;
  W.setName("electrons");
  W.create(1);
  W.createWalkers(1);

  (*W.begin())->Properties(LOCALENERGY) = 1.1;
  (*W.begin())->Properties(LOCALPOTENTIAL) = 1.2;

  std::vector<observable_helper*> h5desc;

  hid_t h_file= H5Fcreate("tmp_obs.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  le_est.registerObservables(h5desc, h_file);
  H5Fclose(h_file);
  // Should make sure h5 file was created?  Check contents?

  LocalEnergyEstimator::RecordListType record;
  le_est.add2Record(record);
  // Not sure how to test this - for now make sure it doesn't crash

}

}
