//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Estimators/EstimatorManagerBase.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"


#include <stdio.h>
#include <string>
#include <random>


using std::string;

namespace qmcplusplus
{
TEST_CASE("Fixed node branch", "[drivers][walker_control]")
{
  Communicate* c = OHMMS::Controller;

  auto emb_uptr = std::make_unique<EstimatorManagerBase>(c);
  auto emb      = emb_uptr.get();

  double tau = 0.5;
  int nideal = 1;
  SimpleFixedNodeBranch fnb(tau, nideal);

  fnb.setEstimatorManager(std::move(emb_uptr));
  REQUIRE(fnb.getEstimatorManager() == emb);

  CHECK(fnb.getTau() == Approx(tau));

  fnb.advanceQMCCounter();
  REQUIRE(fnb.iParam[SimpleFixedNodeBranch::B_COUNTER] == 0);

  // default is classic cutoff scheme
  //fnb.setBranchCutoff();
}

} // namespace qmcplusplus
