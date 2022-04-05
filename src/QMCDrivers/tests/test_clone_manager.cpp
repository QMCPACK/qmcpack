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


#include "Configuration.h"
#include "Message/Communicate.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/QMCUpdateBase.h"
#include "QMCHamiltonians/HamiltonianPool.h"


#include <stdio.h>
#include <string>
#include <random>


using std::string;

namespace qmcplusplus
{
class FakeUpdate : public QMCUpdateBase
{
public:
  FakeUpdate(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator& rg)
      : QMCUpdateBase(w, psi, h, rg)
  {}

  void advanceWalker(Walker_t& thisWalker, bool recompute) override {}
};

TEST_CASE("QMCUpdate", "[drivers]")
{
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);
  elec.setName("e");
  elec.create({1});
  elec.createWalkers(1);

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(massIdx, upIdx)     = 1.0;

  FakeRandom rg;

  QMCHamiltonian h;
  TrialWaveFunction psi;
  FakeUpdate update(elec, psi, h, rg);

  update.put(NULL);

  //update.resetRun(brancher, estimator_manager);
}


TEST_CASE("CloneManager", "[drivers]")
{
  Communicate* c = OHMMS::Controller;

  CloneManager cm;
}

} // namespace qmcplusplus
