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

#include "Configuration.h"
#include "Message/Communicate.h"
#include "QMCDrivers/Crowd.h"
#include "Estimators/tests/FakeEstimator.h"

namespace qmcplusplus
{
TEST_CASE("Crowd integration", "[drivers]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* comm = OHMMS::Controller;

  EstimatorManagerBase em(comm);

  FakeEstimator* fake_est = new FakeEstimator;

  em.add(fake_est, "fake");

  ScalarEstimatorBase* est2 = em.getEstimator("fake");
  FakeEstimator* fake_est2  = dynamic_cast<FakeEstimator*>(est2);
  REQUIRE(fake_est2 != NULL);
  REQUIRE(fake_est2 == fake_est);

  Crowd crowd(em);

  
}
  
}
