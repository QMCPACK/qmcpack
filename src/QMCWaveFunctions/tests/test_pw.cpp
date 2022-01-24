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
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("PlaneWave SPO from HDF for BCC H", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  // BCC H
  PtclOnLatticeTraits::ParticleLayout lattice;
  lattice.R(0, 0) = 3.77945227;
  lattice.R(0, 1) = 0.0;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.77945227;
  lattice.R(1, 2) = 0.0;
  lattice.R(2, 0) = 0.0;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.77945227;
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 1.88972614;
  ions.R[1][1] = 1.88972614;
  ions.R[1][2] = 1.88972614;


  std::vector<int> agroup(2);
  agroup[0] = 1;
  agroup[1] = 1;
  elec.create(agroup);

  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;


  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec.addTable(ions);
  elec.resetGroups();
  elec.update();

  //BCC H
  const char* particles = "<tmp> \
<determinantset type=\"PW\" href=\"bccH.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\"> \
   <slaterdeterminant> \
     <determinant id=\"updet\" size=\"1\"> \
      <occupation mode=\"ground\" spindataset=\"0\"/> \
     </determinant> \
     <determinant id=\"downdet\" size=\"1\"> \
        <occupation mode=\"ground\" spindataset=\"0\"/> \
     </determinant> \
  </slaterdeterminant> \
</determinantset> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr pw1 = xmlFirstElementChild(root);


  PWOrbitalBuilder pw_builder(c, elec, ptcl.getPool());
  auto orb      = pw_builder.buildComponent(pw1);
  SlaterDet* sd = dynamic_cast<SlaterDet*>(orb.get());
  REQUIRE(sd != nullptr);
  REQUIRE(sd->Dets.size() == 2);
  SPOSetPtr spo = sd->getPhi(0);
  REQUIRE(spo != nullptr);
  //SPOSet *spo = einSet.createSPOSetFromXML(ein1);
  //REQUIRE(spo != nullptr);

  int orbSize = spo->getOrbitalSetSize();
  elec.update();
  SPOSet::ValueVector orbs(orbSize);
  spo->evaluateValue(elec, 0, orbs);

  REQUIRE(std::real(orbs[0]) == Approx(-1.2473558998));

#if 0
  // Dump values of the orbitals
  int basisSize= spo->getBasisSetSize();
  printf("orb size = %d basis set size = %d\n",orbSize, basisSize);

  elec.R[1][1] = 0.0;
  double step = 3.78/10;
  FILE *fspo = fopen("spo.dat", "w");
  for (int ix = 0; ix < 10; ix++) {
    for (int iy = 0; iy < 10; iy++) {
      for (int iz = 0; iz < 10; iz++) {
        double x = step*ix;
        double y = step*iy;
        double z = step*iz;
        elec.R[0][0] = x;
        elec.R[0][1] = y;
        elec.R[0][2] = z;
        elec.update();
        SPOSet::ValueVector orbs(orbSize);
        spo->evaluateValue(elec, 0, orbs);
        fprintf(fspo, "%g %g %g",x,y,z);
        for (int j = 0; j < orbSize; j++) {
#ifdef QMC_COMPLEX
          fprintf(fspo, " %g,%g ",orbs[j].real(),orbs[j].imag());
#else
          fprintf(fspo, " %g ",orbs[j]);
#endif
        }
        fprintf(fspo, "\n");
      }
    }
  }
  fclose(fspo);
#endif
}


TEST_CASE("PlaneWave SPO from HDF for LiH arb", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  // LiH
  PtclOnLatticeTraits::ParticleLayout lattice;
  lattice.R(0, 0) = -3.55;
  lattice.R(0, 1) = 0.0;
  lattice.R(0, 2) = 3.55;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.55;
  lattice.R(1, 2) = 3.55;
  lattice.R(2, 0) = -3.55;
  lattice.R(2, 1) = 3.55;
  lattice.R(2, 2) = 0.0;
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 3.55;
  ions.R[1][1] = 3.55;
  ions.R[1][2] = 3.55;


  std::vector<int> agroup(2);
  agroup[0] = 2;
  agroup[1] = 2;
  elec.create(agroup);

  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.R[0][0] = 0.0;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;
  elec.R[2][0] = 0.0;
  elec.R[2][1] = 0.0;
  elec.R[2][2] = 1.0;
  elec.R[3][0] = 0.0;
  elec.R[3][1] = 1.0;
  elec.R[3][2] = 1.0;


  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec.addTable(ions);
  elec.resetGroups();
  elec.update();

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"PW\" href=\"LiH-arb.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\"> \
   <slaterdeterminant> \
     <determinant id=\"updet\" size=\"2\"> \
      <occupation mode=\"ground\" spindataset=\"0\"/> \
     </determinant> \
     <determinant id=\"downdet\" size=\"2\"> \
        <occupation mode=\"ground\" spindataset=\"0\"/> \
     </determinant> \
  </slaterdeterminant> \
</determinantset> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr pw1 = xmlFirstElementChild(root);


  PWOrbitalBuilder pw_builder(c, elec, ptcl.getPool());
  auto orb      = pw_builder.buildComponent(pw1);
  SlaterDet* sd = dynamic_cast<SlaterDet*>(orb.get());
  REQUIRE(sd != nullptr);
  REQUIRE(sd->Dets.size() == 2);
  SPOSetPtr spo = sd->getPhi(0);
  REQUIRE(spo != nullptr);
  //SPOSet *spo = einSet.createSPOSetFromXML(ein1);
  //REQUIRE(spo != nullptr);

  int orbSize = spo->getOrbitalSetSize();
  elec.update();
  SPOSet::ValueVector orbs(orbSize);
  spo->evaluateValue(elec, 0, orbs);

  REQUIRE(std::real(orbs[0]) == Approx(-14.3744302974));

#if 0
  // Dump values of the orbitals
  int basisSize= spo->getBasisSetSize();
  printf("orb size = %d basis set size = %d\n",orbSize, basisSize);

  elec.R[1][1] = 0.0;
  double step = 3.55/10;
  FILE *fspo = fopen("spo.dat", "w");
  for (int ix = 0; ix < 10; ix++) {
    for (int iy = 0; iy < 10; iy++) {
      for (int iz = 0; iz < 10; iz++) {
        double x = step*ix;
        double y = step*iy;
        double z = step*iz;
        elec.R[0][0] = x;
        elec.R[0][1] = y;
        elec.R[0][2] = z;
        elec.update();
        SPOSet::ValueVector orbs(orbSize);
        spo->evaluateValue(elec, 0, orbs);
        fprintf(fspo, "%g %g %g",x,y,z);
        for (int j = 0; j < orbSize; j++) {
#ifdef QMC_COMPLEX
          fprintf(fspo, " %g,%g ",orbs[j].real(),orbs[j].imag());
#else
          fprintf(fspo, " %g ",orbs[j]);
#endif
        }
        fprintf(fspo, "\n");
      }
    }
  }
  fclose(fspo);
#endif
}
} // namespace qmcplusplus
