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


#include "Message/catch_mpi_main.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#ifndef ENABLE_SOA
#include "Particle/SymmetricDistanceTableData.h"
#endif
#include "QMCApp/ParticleSetPool.h"
#include "QMCHamiltonians/BareKineticEnergy.h"

#include "QMCWaveFunctions/TrialWaveFunction.h"
#ifndef ENABLE_SOA
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#endif
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Bare Kinetic Energy", "[hamiltonian]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  elec.create(2);
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 1.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 1.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  int downIdx          = tspecies.addSpecies("d");
  //int chargeIdx = tspecies.addAttribute("charge");
  int massIdx = tspecies.addAttribute("mass");
  //tspecies(chargeIdx, upIdx) = -1;
  //tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

#ifdef ENABLE_SOA
  elec.addTable(ions, DT_SOA);
#else
  elec.addTable(ions, DT_AOS);
#endif
  elec.update();


  ParticleSetPool ptcl = ParticleSetPool(c);


  const char* particles = "<tmp> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr h1 = xmlFirstElementChild(root);

  // This constructor has compile errors (Ps member not initialized)
  //BareKineticEnergy<double> bare_ke;

  BareKineticEnergy<double> bare_ke(elec);
  bare_ke.put(h1);

  elec.L[0] = 1.0;
  elec.L[1] = 0.0;
  double v  = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);

  elec.L[0]    = 0.0;
  elec.L[1]    = 0.0;
  elec.G[0][0] = 1.0;
  elec.G[0][1] = 0.0;
  elec.G[0][2] = 0.0;
  elec.G[1][0] = 0.0;
  elec.G[1][1] = 0.0;
  elec.G[1][2] = 0.0;
  v            = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);
}

TEST_CASE("Bare KE Pulay PBC", "[hamiltonian]")
{
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::PosType PosType;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(20);
  Lattice.LR_dim_cutoff = 15;
  Lattice.reset();


  ParticleSet ions;
  ParticleSet elec;

  ions.setName("ion0");
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 6.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;


  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("Na");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(iatnumber, pIdx)  = 11;
  ions.Lattice = Lattice;
  ions.createSK();


  elec.Lattice = Lattice;
  elec.setName("e");
  elec.create(2);
  elec.R[0][0] = 2.0;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 3.0;
  elec.R[1][1] = 0.0;
  elec.R[1][2] = 0.0;


  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

  elec.createSK();

  ParticleSetPool ptcl = ParticleSetPool(c);

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  //Cool.  Now to construct a wavefunction with 1 and 2 body jastrow (no determinant)
  TrialWaveFunction psi(c);

  //Add the two body jastrow
  const char* particles = "<tmp> \
  <jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\">  \
      <correlation speciesA=\"u\" speciesB=\"d\" rcut=\"10\" size=\"8\"> \
          <coefficients id=\"ud\" type=\"Array\"> 2.015599059 1.548994099 1.17959447 0.8769687661 0.6245736507 0.4133517767 0.2333851935 0.1035636904</coefficients> \
        </correlation> \
  </jastrow> \
  </tmp> \
  ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas2 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow(c, elec);
  psi.addComponent(jastrow.buildComponent(jas2), "RadialJastrow");
  // Done with two body jastrow.

  //Add the one body jastrow.
  const char* particles2 = "<tmp> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" source=\"ion0\" print=\"yes\"> \
        <correlation elementType=\"Na\" rcut=\"10\" size=\"10\" cusp=\"0\"> \
          <coefficients id=\"eNa\" type=\"Array\"> 1.244201343 -1.188935609 -1.840397253 -1.803849126 -1.612058635 -1.35993202 -1.083353212 -0.8066295188 -0.5319252448 -0.3158819772</coefficients> \
        </correlation> \
      </jastrow> \
  </tmp> \
  ";
  bool okay3             = doc.parseFromString(particles2);
  REQUIRE(okay3);

  root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow1bdy(c, elec, ions);
  psi.addComponent(jastrow1bdy.buildComponent(jas1), "RadialJastrow");

  const char* kexml = "<tmp> \
</tmp> \
";

  root = doc.getRoot();

  xmlNodePtr h1 = xmlFirstElementChild(root);

  BareKineticEnergy<double> bare_ke(elec);
  bare_ke.put(h1);

  // update all distance tables
  ions.update();
  elec.update();

  RealType logpsi = psi.evaluateLog(elec);

  RealType keval = bare_ke.evaluate(elec);

  //This is validated against an alternate code path (waveefunction tester for local energy).
  REQUIRE(keval == Approx(-0.147507745));

#ifdef ENABLE_SOA
  ParticleSet::ParticlePos_t HFTerm, PulayTerm;
  HFTerm.resize(ions.getTotalNum());
  PulayTerm.resize(ions.getTotalNum());

  RealType keval2 = bare_ke.evaluateWithIonDerivs(elec, ions, psi, HFTerm, PulayTerm);

  REQUIRE(keval2 == Approx(-0.147507745));
  //These are validated against finite differences (delta=1e-6).
  REQUIRE(PulayTerm[0][0] == Approx(-0.13166));
  REQUIRE(PulayTerm[0][1] == Approx(0.0));
  REQUIRE(PulayTerm[0][2] == Approx(0.0));
  REQUIRE(PulayTerm[1][0] == Approx(-0.12145));
  REQUIRE(PulayTerm[1][1] == Approx(0.0));
  REQUIRE(PulayTerm[1][2] == Approx(0.0));

#endif
}
} // namespace qmcplusplus
