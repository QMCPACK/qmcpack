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
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Particle/SymmetricDistanceTableData.h"
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
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  EstimatorManagerBase emb(c);

  double tau = 0.5;
  int nideal = 1;
  SimpleFixedNodeBranch fnb(tau, nideal);

  fnb.setEstimatorManager(&emb);
  REQUIRE(fnb.getEstimatorManager() == &emb);

  REQUIRE(fnb.getTau() == Approx(tau));

  fnb.advanceQMCCounter();
  REQUIRE(fnb.iParam[SimpleFixedNodeBranch::B_COUNTER] == 0);

  // default is classic cutoff scheme
  //fnb.setBranchCutoff();
}

}
