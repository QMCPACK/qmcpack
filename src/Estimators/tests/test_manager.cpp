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
#include "Estimators/EstimatorManagerBase.h"
#include "Estimators/ScalarEstimatorBase.h"


#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

class FakeEstimator : public ScalarEstimatorBase
{
  virtual void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last , RealType wgt) {}

  virtual void add2Record(RecordNamedProperty<RealType>& record) {}

  virtual void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid) {}

  virtual ScalarEstimatorBase* clone() { return new FakeEstimator; }
};

TEST_CASE("EstimatorManagerBase", "[estimators]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  EstimatorManagerBase em(c);

  REQUIRE(em.size() == 0);

  // Must create on heap since the EstimatorManager destructor deletes all estimators
  FakeEstimator *fake_est = new FakeEstimator;

  em.add(fake_est, "fake");

  ScalarEstimatorBase *est2 = em.getEstimator("fake");
  FakeEstimator *fake_est2 = dynamic_cast<FakeEstimator *>(est2);
  REQUIRE(fake_est2 != NULL);
  REQUIRE(fake_est2 == fake_est);

  // Check the copy constructor
  EstimatorManagerBase em2(em);
  REQUIRE(em.size() == 1);

  em.start(2, true);

  em.stop();
  // compute averages over threads
  //em.stop();
  em.reset();

}

}
