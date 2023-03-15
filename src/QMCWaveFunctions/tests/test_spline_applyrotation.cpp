//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//
// File created by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2C.h"
//#include "Utilities/for_testing/checkMatrix.hpp"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{

TEST_CASE("Spline applyRotation no rotation", "[wavefunction]")
{
  // How to get rid of all this annoying boilerplate?
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // monoO
  // lattice.R(0,0) = {5.10509515, -3.23993545,  0.0, 5.10509515, 3.23993545, 0.0, -6.49690625, 0.0, 7.08268015};

  // diamondC_1x1x1
  lattice.R = {3.37316115, 3.37316115, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // Load diamondC_1x1x1 wfn and explicitly construct a SplineC2C object with 7 orbitals
  // This results in padding of the spline coefs table and thus is a more stringent test.
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  /*
      Here we test the low-level spline implementation of applyRotation().
      Manually encode the unitary transformation. Ugly, but it works.

      NB: This is truncated to 5 sig-figs, so there is some slop here as
      compared to what is done in the splines via apply_rotation().
      So below we reduce the threshold for comparison. This can
      probably be ditched once we have a way to grab the actual
      rotation matrix...
    */
  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // 0.) Check the bare splines
  // value
  // TODO: Use checkMatrix(phimat, spomat);

  CHECK(std::real(psiM_bare[1][0]) == Approx(-0.8886948824));
  CHECK(std::real(psiM_bare[1][1]) == Approx(1.4194120169));
  // grad
  CHECK(std::real(dpsiM_bare[1][0][0]) == Approx(-0.0000183403));
  CHECK(std::real(dpsiM_bare[1][0][1]) == Approx(0.1655139178));
  CHECK(std::real(dpsiM_bare[1][0][2]) == Approx(-0.0000193077));
  CHECK(std::real(dpsiM_bare[1][1][0]) == Approx(-1.3131694794));
  CHECK(std::real(dpsiM_bare[1][1][1]) == Approx(-1.1174004078));
  CHECK(std::real(dpsiM_bare[1][1][2]) == Approx(-0.8462534547));
  // lapl
  CHECK(std::real(d2psiM_bare[1][0]) == Approx(1.3313053846));
  CHECK(std::real(d2psiM_bare[1][1]) == Approx(-4.712583065));

  /*
      2.) Apply a rotation to the splines. In this case we apply the identity matrix = no rotation.
    */
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat       = 0;
  rot_mat[0][0] = 1.;
  rot_mat[1][1] = 1.;
  rot_mat[2][2] = 1.;
  rot_mat[3][3] = 1.;
  rot_mat[4][4] = 1.;
  rot_mat[5][5] = 1.;
  rot_mat[6][6] = 1.;

  spo->storeParamsBeforeRotation(); // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat, false);

  // 3.) Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // 4.) Compute the expected values by hand using the transformation above
  SPOSet::ValueType val1{0.};
  SPOSet::ValueType val2{0.};
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    val1 += psiM_bare[0][i] * rot_mat[i][0];
    val2 += psiM_bare[1][i] * rot_mat[i][0];
  }

  // 5.) check that the rotation was done correctly
  // value
  CHECK(std::real(psiM_rot[0][0]) == Approx(std::real(val1)));
  CHECK(std::imag(psiM_rot[0][0]) == Approx(std::imag(val1)));
  CHECK(std::real(psiM_rot[1][0]) == Approx(std::real(val2)));
  CHECK(std::imag(psiM_rot[1][0]) == Approx(std::imag(val2)));

  SPOSet::GradType grad1{0.};
  SPOSet::GradType grad2{0.};
  for (auto j = 0; j < grad1.size(); j++)
  {
    for (auto i = 0; i < rot_mat.size1(); i++)
    {
      grad1[j] += dpsiM_bare[0][i][j] * rot_mat[i][0];
      grad2[j] += dpsiM_bare[1][i][j] * rot_mat[i][0];
    }
  }

  // grad
  CHECK(std::real(dpsiM_rot[0][0][0]) == Approx(std::real(grad1[0])).epsilon(0.0001));
  CHECK(std::real(dpsiM_rot[0][0][1]) == Approx(std::real(grad1[1])).epsilon(0.0001));
  CHECK(std::real(dpsiM_rot[0][0][2]) == Approx(std::real(grad1[2])).epsilon(0.0001));
  CHECK(std::real(dpsiM_rot[1][0][0]) == Approx(std::real(grad2[0])).epsilon(0.0001));
  CHECK(std::real(dpsiM_rot[1][0][1]) == Approx(std::real(grad2[1])).epsilon(0.0001));
  CHECK(std::real(dpsiM_rot[1][0][2]) == Approx(std::real(grad2[2])).epsilon(0.0001));

  SPOSet::ValueType lap1 = 0.;
  SPOSet::ValueType lap2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    lap1 += d2psiM_bare[0][i] * rot_mat[i][0];
    lap2 += d2psiM_bare[1][i] * rot_mat[i][0];
  }

  // Lapl
  CHECK(std::real(d2psiM_rot[0][0]) == Approx(std::real(lap1)).epsilon(0.0001));
  CHECK(std::real(d2psiM_rot[1][0]) == Approx(std::real(lap2)).epsilon(0.0001));

} // TEST_CASE


} // namespace qmcplusplus
