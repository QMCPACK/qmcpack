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
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "PlaneWave/PWOrbitalSetBuilder.h"


#include <stdio.h>
#include <cmath>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("PlaneWave SPO from HDF for BCC H", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  // BCC H
  Lattice lattice;
  lattice.R = {3.77945227, 0.0, 0.0, 0.0, 3.77945227, 0.0, 0.0, 0.0, 3.77945227};
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.createSimulationCellByLattice(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create({2});
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.R[1] = {1.88972614, 1.88972614, 1.88972614};

  elec.create({1, 1});
  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.0, 1.0, 0.0};

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  const int upIdx              = tspecies.addSpecies("u");
  const int downIdx            = tspecies.addSpecies("d");
  const int chargeIdx          = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec.addTable(ions);
  elec.resetGroups();
  elec.update();

  //BCC H
  const char* particles = R"(
<sposet_collection type="PW" href="bccH.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion">
  <sposet name="updet" size="1" spindataset="0">
    <occupation mode="ground"/>
  </sposet>
</sposet_collection>
)";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(particles));

  xmlNodePtr root = doc.getRoot();
  xmlNodePtr pw1  = xmlFirstElementChild(root);


  PWOrbitalSetBuilder pw_builder(elec, c, root);
  auto spo = pw_builder.createSPOSet(pw1);
  REQUIRE(spo);

  const int orbSize = spo->getOrbitalSetSize();
  elec.update();
  SPOSet::ValueVector orbs(orbSize);
  spo->evaluateValue(elec, 0, orbs);

  CHECK(std::real(orbs[0]) == Approx(-1.2473558998));

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
        elec.R[0] = {x, y, z};
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

TEST_CASE("PlaneWave tiled SPO from HDF for diamond C", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  Lattice lattice;
  lattice.R = {6.74632230, 6.74632230, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};
  lattice.reset();

  ParticleSetPool ptcl(c);
  ptcl.createSimulationCellByLattice(lattice);
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& elec(*elec_uptr);
  elec.create({8, 8});
  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.R[0] = {0.25, 0.5, 0.75};
  elec.update();

  const char* particles = R"(
<sposet_collection type="PW" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0">
  <sposet name="updet" size="8" spindataset="0"/>
</sposet_collection>
)";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(particles));

  xmlNodePtr root = doc.getRoot();
  PWOrbitalSetBuilder pw_builder(elec, c, root);
  auto spo = pw_builder.createSPOSet(xmlFirstElementChild(root));
  REQUIRE(spo);
  CHECK(spo->getClassName() == "CompositeSPOSet");
  REQUIRE(spo->getOrbitalSetSize() == 8);

  SPOSet::ValueVector values(8);
  SPOSet::GradVector gradients(8);
  SPOSet::ValueVector laplacians(8);
  spo->evaluateVGL(elec, 0, values, gradients, laplacians);

  SPOSet::ValueMatrix value_matrix(1, 8);
  SPOSet::GradMatrix gradient_matrix(1, 8);
  SPOSet::HessMatrix hessian_matrix(1, 8);
  spo->evaluate_notranspose(elec, 0, 1, value_matrix, gradient_matrix, hessian_matrix);
#if defined(MIXED_PRECISION)
  constexpr double hessian_tolerance = 1e-4;
#else
  constexpr double hessian_tolerance = 1e-10;
#endif
  for (int orbital_index = 0; orbital_index < 8; ++orbital_index)
  {
    CHECK(std::isfinite(std::abs(values[orbital_index])));
    CHECK(std::isfinite(std::abs(laplacians[orbital_index])));
    CHECK(std::isfinite(std::abs(gradients[orbital_index][0])));
    CHECK(std::abs(value_matrix(0, orbital_index) - values[orbital_index]) == Approx(0.0).margin(1e-12));
    SPOSet::ValueType hessian_trace = 0;
    for (int row = 0; row < OHMMS_DIM; ++row)
    {
      hessian_trace += hessian_matrix(0, orbital_index)(row, row);
      for (int column = 0; column < OHMMS_DIM; ++column)
        CHECK(std::isfinite(std::abs(hessian_matrix(0, orbital_index)(row, column))));
    }
    CHECK(std::abs(hessian_trace - laplacians[orbital_index]) == Approx(0.0).margin(hessian_tolerance));
  }
}


TEST_CASE("PlaneWave SPO from HDF for LiH arb", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  // LiH
  Lattice lattice;
  lattice.R = {-3.55, 0.0, 3.55, 0.0, 3.55, 3.55, -3.55, 3.55, 0.0};
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.createSimulationCellByLattice(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create({1, 1});
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.R[1] = {3.55, 3.55, 3.55};

  elec.create({2, 2});
  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.0, 1.0, 0.0};
  elec.R[2] = {0.0, 0.0, 1.0};
  elec.R[3] = {0.0, 1.0, 1.0};

  SpeciesSet& tspecies         = elec.getSpeciesSet();
  const int upIdx              = tspecies.addSpecies("u");
  const int downIdx            = tspecies.addSpecies("d");
  const int chargeIdx          = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec.addTable(ions);
  elec.resetGroups();
  elec.update();

  //diamondC_1x1x1
  const char* particles = R"(
<sposet_collection type="PW" href="LiH-arb.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion">
  <sposet name="updet" size="2" spindataset="0">
    <occupation mode="ground"/>
  </sposet>
</sposet_collection>
)";

  Libxml2Document doc;
  REQUIRE(doc.parseFromString(particles));

  xmlNodePtr root = doc.getRoot();
  xmlNodePtr pw1  = xmlFirstElementChild(root);


  PWOrbitalSetBuilder pw_builder(elec, c, root);
  auto spo = pw_builder.createSPOSet(pw1);
  REQUIRE(spo);

  const int orbSize = spo->getOrbitalSetSize();
  elec.update();
  SPOSet::ValueVector orbs(orbSize);
  spo->evaluateValue(elec, 0, orbs);

  CHECK(std::real(orbs[0]) == Approx(-14.3744302974));

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
        elec.R[0] = {x, y, z};
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
