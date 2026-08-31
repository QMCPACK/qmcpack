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
#include "PlaneWave/PWBandInfo.h"
#include "PlaneWave/PWOrbitalSetBuilder.h"
#include "PlaneWave/PWTiling.h"


#include <stdio.h>
#include <cmath>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("PlaneWave band sort modes", "[wavefunction][plane_wave]")
{
  const std::vector<pw::BandInfo> input{{2.0, 0, 0}, {1.0, 0, 1}, {2.0 + 0.5e-6, 1, 0}, {1.0 + 0.5e-6, 1, 1}};

  auto band_indices = [](const std::vector<pw::BandInfo>& bands) {
    std::vector<std::pair<int, int>> indices;
    for (const pw::BandInfo& band : bands)
      indices.emplace_back(band.twist_index, band.band_index);
    return indices;
  };

  SECTION("insertion order")
  {
    auto bands = input;
    pw::sortBands(bands, 0);
    CHECK(band_indices(bands) == std::vector<std::pair<int, int>>{{0, 0}, {0, 1}, {1, 0}, {1, 1}});
  }
  SECTION("energy order")
  {
    auto bands = input;
    pw::sortBands(bands, 1);
    CHECK(band_indices(bands) == std::vector<std::pair<int, int>>{{0, 1}, {1, 1}, {0, 0}, {1, 0}});
  }
  SECTION("band index order")
  {
    auto bands = input;
    pw::sortBands(bands, 2);
    CHECK(band_indices(bands) == std::vector<std::pair<int, int>>{{0, 0}, {1, 0}, {0, 1}, {1, 1}});
  }
  SECTION("invalid order")
  {
    auto bands = input;
    CHECK_THROWS_AS(pw::sortBands(bands, 3), std::runtime_error);
  }
}

TEST_CASE("PlaneWave tilematrix folding", "[wavefunction][plane_wave]")
{
  pw::TileMatrix tile_matrix;
  for (int row = 0; row < OHMMS_DIM; ++row)
    for (int column = 0; column < OHMMS_DIM; ++column)
      tile_matrix(row, column) = 0;
  tile_matrix(0, 0) = 1;
  tile_matrix(0, 1) = 1;
  tile_matrix(1, 0) = 1;
  tile_matrix(1, 1) = -1;
  tile_matrix(2, 2) = 1;

  const std::vector<pw::Twist> primitive_twists{{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}};
  const pw::TwistGroups groups = pw::groupPrimitiveTwists(tile_matrix, primitive_twists);
  REQUIRE(groups.super_twists.size() == 1);
  CHECK(groups.primitive_indices[0] == std::vector<int>{0, 1});
  CHECK_NOTHROW(pw::validateTwistGroup(tile_matrix, primitive_twists, groups, 0));
  CHECK(pw::findTwistGroup(groups, {1.0, 0.0, 0.0}) == 0);
  CHECK(pw::canonicalizeTwist({0.6, 0.0, 0.0})[0] == Approx(-0.4));

  const pw::Lattice primitive_lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  const pw::Lattice super_lattice = pw::makeSuperLattice(tile_matrix, primitive_lattice);
  CHECK_NOTHROW(pw::validateLattice(tile_matrix, primitive_lattice, super_lattice));
  pw::Lattice mismatched_lattice = super_lattice;
  mismatched_lattice(0, 0) += 0.1;
  CHECK_THROWS_AS(pw::validateLattice(tile_matrix, primitive_lattice, mismatched_lattice), std::runtime_error);

  const std::vector<pw::Twist> incomplete_twists{{0.0, 0.0, 0.0}};
  const pw::TwistGroups incomplete_groups = pw::groupPrimitiveTwists(tile_matrix, incomplete_twists);
  CHECK_THROWS_AS(pw::validateTwistGroup(tile_matrix, incomplete_twists, incomplete_groups, 0), std::runtime_error);

  const std::vector<pw::Twist> duplicate_twists{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  const pw::TwistGroups duplicate_groups = pw::groupPrimitiveTwists(tile_matrix, duplicate_twists);
  CHECK_THROWS_AS(pw::validateTwistGroup(tile_matrix, duplicate_twists, duplicate_groups, 0), std::runtime_error);

  pw::TileMatrix singular_matrix;
  for (int row = 0; row < OHMMS_DIM; ++row)
    for (int column = 0; column < OHMMS_DIM; ++column)
      singular_matrix(row, column) = 0;
  CHECK_THROWS_AS(pw::validateTwistGroup(singular_matrix, incomplete_twists, incomplete_groups, 0), std::runtime_error);
}

TEST_CASE("PlaneWave real twist pairing", "[wavefunction][plane_wave]")
{
  const std::vector<pw::Twist> twists{{-1.0 / 3.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0 / 3.0, 0.0, 0.0}};
  const std::vector<pw::DistinctTwist> distinct = pw::findDistinctRealTwists({0, 1, 2}, twists);
  REQUIRE(distinct.size() == 2);
  CHECK(distinct[0].twist_index == 1);
  CHECK_FALSE(distinct[0].make_two_copies);
  CHECK(distinct[1].twist_index == 2);
  CHECK(distinct[1].make_two_copies);
  CHECK(pw::isRealCompatible({0.5, 0.0, 0.0}));
  CHECK_FALSE(pw::isRealCompatible({1.0 / 3.0, 0.0, 0.0}));
}

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
  CHECK(dot(elec.getTwist(), elec.getTwist()) == Approx(0.0));

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

  const char* smaller_spo_input = R"(
<sposet_collection type="PW" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0">
  <sposet name="small_spo" size="3" spindataset="0"/>
</sposet_collection>
)";
  Libxml2Document smaller_doc;
  REQUIRE(smaller_doc.parseFromString(smaller_spo_input));
  xmlNodePtr smaller_root = smaller_doc.getRoot();
  PWOrbitalSetBuilder smaller_builder(elec, c, smaller_root);
  auto smaller_spo = smaller_builder.createSPOSet(xmlFirstElementChild(smaller_root));
  REQUIRE(smaller_spo);
  CHECK(smaller_spo->getOrbitalSetSize() == 3);
}


TEST_CASE("PlaneWave SPO from HDF for LiH X", "[wavefunction]")
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

  const char* particles = R"(
<sposet_collection type="PW" href="LiH-x.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="-1" twist="0.5 0 0" source="ion">
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
  CHECK(elec.getTwist()[0] == Approx(0.5));
  CHECK(elec.getTwist()[1] == Approx(0.0));
  CHECK(elec.getTwist()[2] == Approx(0.0));

  const int orbSize = spo->getOrbitalSetSize();
  elec.update();
  SPOSet::ValueVector orbs(orbSize);
  spo->evaluateValue(elec, 0, orbs);

  CHECK(std::isfinite(std::abs(orbs[0])));

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
