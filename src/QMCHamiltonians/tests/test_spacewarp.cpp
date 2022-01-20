//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "QMCHamiltonians/SpaceWarpTransformation.h"
//#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
//#include "QMCHamiltonians/ForceChiesaPBCAA.h"
//#include "QMCHamiltonians/ForceCeperley.h"
//#include "QMCHamiltonians/CoulombPotential.h"
//#include "QMCHamiltonians/CoulombPBCAA.h"
//#include "QMCHamiltonians/CoulombPBCAB.h"
//#include "QMCWaveFunctions/TrialWaveFunction.h"
//#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
//#include "QMCWaveFunctions/Fermion/SlaterDet.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using RealType = QMCTraits::RealType;
TEST_CASE("SpaceWarp", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("Na2.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  const SimulationCell simulation_cell;

  ParticleSet ions(simulation_cell);
  XMLParticleParser parse_ions(ions);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.put(particleset_ion[0]);

  REQUIRE(ions.groups() == 1);
  REQUIRE(ions.R.size() == 2);
  ions.update();

  ParticleSet elec(simulation_cell);
  XMLParticleParser parse_elec(elec);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.put(particleset_elec[0]);

  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 2);

  elec.addTable(ions);
  elec.update();

  //Now build the wavefunction.  This will be needed to test \Nabla_i E_L and \Nabla_i logPsi contributions.
  //For now, just take them from a reference calculation.

  using Force_t = ParticleSet::ParticlePos_t;
  Force_t dE_L;
  Force_t el_contribution;
  Force_t psi_contribution;

  dE_L.resize(elec.getTotalNum());
  el_contribution.resize(ions.getTotalNum());
  psi_contribution.resize(ions.getTotalNum());

  dE_L[0][0] = -0.0328339806050;
  dE_L[0][1] = -0.0834441565340;
  dE_L[0][2] = 0.0997813066140;
  dE_L[1][0] = -0.0140597469190;
  dE_L[1][1] = 0.0591827022730;
  dE_L[1][2] = -0.0622852142310;

  elec.G[0][0] = 0.4167938814700;
  elec.G[0][1] = 0.2878426639600;
  elec.G[0][2] = -0.3470187402100;
  elec.G[1][0] = -0.2946265813200;
  elec.G[1][1] = -0.3606166249000;
  elec.G[1][2] = 0.3751881159300;

  SpaceWarpTransformation swt(elec, ions);
  swt.setPow(3.0);

  REQUIRE(swt.f(2.0) == Approx(0.125));
  REQUIRE(swt.df(2.0) == Approx(-0.1875));

  swt.setPow(4.0);
  REQUIRE(swt.f(2.0) == Approx(0.0625));
  REQUIRE(swt.df(2.0) == Approx(-0.125));

  swt.computeSWT(elec, ions, dE_L, elec.G, el_contribution, psi_contribution);
  app_log() << "EL_Contribution:  " << el_contribution << std::endl;
  app_log() << "PSi_Contribution: " << psi_contribution << std::endl;
  REQUIRE(el_contribution[0][0] == Approx(-0.0326934696861));
  REQUIRE(el_contribution[0][1] == Approx(-0.0826080664130));
  REQUIRE(el_contribution[0][2] == Approx(0.0988243408507));
  REQUIRE(el_contribution[1][0] == Approx(-0.0142002578379));
  REQUIRE(el_contribution[1][1] == Approx(0.0583466121520));
  REQUIRE(el_contribution[1][2] == Approx(-0.0613282484677));

  REQUIRE(psi_contribution[0][0] == Approx(0.4051467191368));
  REQUIRE(psi_contribution[0][1] == Approx(0.2757724717133));
  REQUIRE(psi_contribution[0][2] == Approx(-0.3334287440127));
  REQUIRE(psi_contribution[1][0] == Approx(-0.2829794189868));
  REQUIRE(psi_contribution[1][1] == Approx(-0.3485464326533));
  REQUIRE(psi_contribution[1][2] == Approx(0.3615981197327));
}
} //namespace qmcplusplus
