//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/DensityEstimator.h"

#include <stdio.h>
#include <string>


namespace qmcplusplus
{
TEST_CASE("Density Estimator", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.setName("elec");
  elec.create({2});
  elec.R[0] = {0.0, 1.0, 0.0};
  elec.R[1] = {0.4, 0.3, 0.0};

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx                 = tspecies.addAttribute("mass");
  int eChargeIdx              = tspecies.addAttribute("charge");
  tspecies(eChargeIdx, upIdx) = -1.0;
  tspecies(massIdx, upIdx)    = 1.0;

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();
  elec.update();

  DensityEstimator density_estimator(elec);

  ParticleSet::Walker_t dummy = ParticleSet::Walker_t(1);
  density_estimator.setHistories(dummy);

  const double val = density_estimator.evaluate(elec);
  CHECK(val == Approx(0.0));
}

TEST_CASE("Density Estimator evaluate exception", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);
  DensityEstimator density_estimator(elec);
  REQUIRE(density_estimator.getClassName() == "DensityEstimator");
  // empty function
  density_estimator.resetTargetParticleSet(elec);
  // check exception is thrown when weights are not defined
  CHECK_THROWS_AS(density_estimator.evaluate(elec), std::invalid_argument);
}

} // namespace qmcplusplus
