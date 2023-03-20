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
#include "Utilities/for_testing/checkMatrix.hpp"

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

  // diamondC_1x1x1
  ParticleSet::ParticleLayout lattice;
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
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueType val{0.};
  for (auto i = 0; i < elec_.R.size(); i++)
    {
      for (auto j = 0; j < orbitalsetsize; j++)
        {
          psiM_rot_manual[i][j] = 0.;
          val = 0.;
          for (auto k = 0; k < orbitalsetsize; k++)
            {
              val += psiM_bare[i][k] * rot_mat[k][j];
            }
          psiM_rot_manual[i][j] = val;
        }
  }

  // 5.) check that the rotation was done correctly
  // value
  std::cerr << "JPT DEBUG :: check value...\n";
  checkMatrix(psiM_rot_manual, psiM_rot);

  // grad
  std::cerr << "JPT DEBUG :: check grad...\n";
  SPOSet::GradMatrix grad_rot_manual(elec_.R.size(), orbitalsetsize);
  std::vector<SPOSet::ValueType> grad_val(3);
  for (auto i = 0; i < elec_.R.size(); i++)
  {
    for (auto j = 0; j < orbitalsetsize; j++)
      {
        grad_rot_manual[i][j][0] = 0.;
        grad_rot_manual[i][j][1] = 0.;
        grad_rot_manual[i][j][2] = 0.;
        std::fill(grad_val.begin(), grad_val.end(), 0);
        for (auto k = 0; k < orbitalsetsize; k++)
          {
            for (auto l = 0; l < 3; l++)
              {
                grad_val[l] += dpsiM_bare[i][k][l] * rot_mat[k][j];
              }
          }
        grad_rot_manual[i][j][0] = grad_val[0];
        grad_rot_manual[i][j][1] = grad_val[1];
        grad_rot_manual[i][j][2] = grad_val[2];
      }
  }

  // Ugh... no checkMatrix for tensors?
  //checkMatrix(grad_rot_manual, dpsiM_rot);
  SPOSet::ValueType res = 0.;
  for (auto i=0; i<elec_.R.size(); i++ )
    {
      for (auto j=0; j<orbitalsetsize; j++ )
        {
          for (auto k=0; k<3; k++ )
            {
              res += std::norm(grad_rot_manual[i][j][k] - dpsiM_rot[i][j][k]);
            }
        }
    }
  CHECK(std::real(res) == Approx(0.));
  CHECK(std::imag(res) == Approx(0.));
  std::cerr << "JPT DEBUG: res = " << res << "\n";

  /*
  // laplacian
  SPOSet::ValueType lap1 = 0.;
  SPOSet::ValueType lap2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    lap1 += d2psiM_bare[0][i] * rot_mat[i][0];
    lap2 += d2psiM_bare[1][i] * rot_mat[i][0];
  }
  CHECK(std::real(d2psiM_rot[0][0]) == Approx(std::real(lap1)).epsilon(0.0001));
  CHECK(std::real(d2psiM_rot[1][0]) == Approx(std::real(lap2)).epsilon(0.0001));
  */
} // TEST_CASE


} // namespace qmcplusplus
