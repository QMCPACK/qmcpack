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
TEST_CASE("Spline applyRotation zero rotation", "[wavefunction]")
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
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" gpu="no" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Check before rotation is applied. From 'test_einset.cpp'
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

  // zero rotation = identity matrix
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat = 0;
  for (int i = 0; i < orbitalsetsize; i++)
    rot_mat[i][i] = 1;

  spo->storeParamsBeforeRotation(); // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat, false);

  // Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Compute the expected value, gradient, and laplacian by hand using the transformation above
  // Check value
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        psiM_rot_manual[i][j] += psiM_bare[i][k] * rot_mat[k][j];
    }
  checkMatrix(psiM_rot_manual, psiM_rot);

  // Check grad
  SPOSet::GradMatrix dpsiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      dpsiM_rot_manual[i][j][0] = 0.;
      dpsiM_rot_manual[i][j][1] = 0.;
      dpsiM_rot_manual[i][j][2] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        for (int l = 0; l < 3; l++)
          dpsiM_rot_manual[i][j][l] += dpsiM_bare[i][k][l] * rot_mat[k][j];
    }

  // No checkMatrix for tensors? Gotta do it by hand...
  double res = 0.;
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
      for (int k = 0; k < 3; k++)
        res += std::sqrt(std::norm(dpsiM_rot_manual[i][j][k] - dpsiM_rot[i][j][k]));

  CHECK(res == Approx(0.).epsilon(2e-4)); // increase threshold to accomodate mixed precision.

  // Check laplacian
  SPOSet::ValueMatrix d2psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      d2psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        d2psiM_rot_manual[i][j] += d2psiM_bare[i][k] * rot_mat[k][j];
    }

  checkMatrix(d2psiM_rot_manual, d2psiM_rot);

} // TEST_CASE


TEST_CASE("Spline applyRotation one rotation", "[wavefunction]")
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
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" gpu="no" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Check before rotation is applied. From 'test_einset.cpp'
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

  /* Apply a rotation, U = exp(-K) as shown below, to the splines.
     Because K = -K.T, U is unitary.
     K=
      0.00000000e+00  -2.33449314e-01   5.21735139e-01   1.67744276e-01  -1.68565994e-01   1.26632312e-02   2.29272040e-01
      2.33449314e-01   0.00000000e+00  -3.81982621e-01  -6.76995896e-01   7.15727948e-02   9.10926961e-01  -6.94864205e-01
     -5.21735139e-01   3.81982621e-01   0.00000000e+00   2.37433139e-01  -7.09878105e-01  -6.33841289e-01  -7.49009582e-01
     -1.67744276e-01   6.76995896e-01  -2.37433139e-01   0.00000000e+00  -6.72002172e-01   4.27152921e-01  -4.82983743e-03
      1.68565994e-01  -7.15727948e-02   7.09878105e-01   6.72002172e-01   0.00000000e+00  -7.93860012e-01  -7.76624731e-01
     -1.26632312e-02  -9.10926961e-01   6.33841289e-01  -4.27152921e-01   7.93860012e-01   0.00000000e+00  -2.74931052e-01
     -2.29272040e-01   6.94864205e-01   7.49009582e-01   4.82983743e-03   7.76624731e-01   2.74931052e-01   0.00000000e+00
     U=
      8.06061880e-01   3.54921598e-01  -3.40706426e-01  -5.90163619e-02   1.11650454e-01  -1.99768450e-01  -2.28818375e-01
      1.58069821e-02   2.10363421e-01   2.74922448e-01   2.12581764e-01   2.64602356e-01  -6.34971914e-01   6.01265560e-01
      4.43696646e-01  -6.06912539e-02   2.61413193e-01  -1.98368802e-01  -6.43234645e-02   5.75880430e-01   5.96646495e-01
      1.45865363e-01  -5.16577220e-01  -2.09866367e-01   5.55699395e-01   5.73062990e-01   1.74778224e-01   8.77170506e-03
     -3.24609748e-01   2.89664179e-01  -7.50613752e-01  -1.63060005e-01   1.70803377e-01   1.63784167e-01   4.05850414e-01
      1.32570771e-03   3.80299846e-01  -5.08439810e-02   7.59141791e-01  -4.77844928e-01   2.12149087e-01   5.60882349e-02
      1.62781208e-01  -5.75073150e-01  -3.60485665e-01  -1.70070331e-02  -5.72251258e-01  -3.50549638e-01   2.49394158e-01
     UU.T=
      1.00000000e+00  -2.77555756e-17   5.55111512e-17  -8.67361738e-18  -1.11022302e-16  -6.07153217e-17   6.93889390e-17
     -2.77555756e-17   1.00000000e+00   5.55111512e-17  -1.84748050e-16   5.55111512e-17  -6.93889390e-18   8.32667268e-17
      5.55111512e-17   5.55111512e-17   1.00000000e+00   1.45716772e-16  -5.55111512e-17   1.38777878e-16   5.55111512e-17
     -8.67361738e-18  -1.84748050e-16   1.45716772e-16   1.00000000e+00  -1.95156391e-17   1.80411242e-16   8.02309608e-17
     -1.11022302e-16   5.55111512e-17  -5.55111512e-17  -1.95156391e-17   1.00000000e+00   7.28583860e-17   1.11022302e-16
     -6.07153217e-17  -6.93889390e-18   1.38777878e-16   1.80411242e-16   7.28583860e-17   1.00000000e+00  -1.12757026e-16
      6.93889390e-17   8.32667268e-17   5.55111512e-17   8.02309608e-17   1.11022302e-16  -1.12757026e-16   1.00000000e+00

      NB: There's nothing special about this rotation. I purposefully made something 'ugly' to make a worst case scenario...
    */
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat[0][0] = 8.06061880e-01;
  rot_mat[0][1] = 3.54921598e-01;
  rot_mat[0][2] = -3.40706426e-01;
  rot_mat[0][3] = -5.90163619e-02;
  rot_mat[0][4] = 1.11650454e-01;
  rot_mat[0][5] = -1.99768450e-01;
  rot_mat[0][6] = -2.28818375e-01;
  rot_mat[1][0] = 1.58069821e-02;
  rot_mat[1][1] = 2.10363421e-01;
  rot_mat[1][2] = 2.74922448e-01;
  rot_mat[1][3] = 2.12581764e-01;
  rot_mat[1][4] = 2.64602356e-01;
  rot_mat[1][5] = -6.34971914e-01;
  rot_mat[1][6] = 6.01265560e-01;
  rot_mat[2][0] = 4.43696646e-01;
  rot_mat[2][1] = -6.06912539e-02;
  rot_mat[2][2] = 2.61413193e-01;
  rot_mat[2][3] = -1.98368802e-01;
  rot_mat[2][4] = -6.43234645e-02;
  rot_mat[2][5] = 5.75880430e-01;
  rot_mat[2][6] = 5.96646495e-01;
  rot_mat[3][0] = 1.45865363e-01;
  rot_mat[3][1] = -5.16577220e-01;
  rot_mat[3][2] = -2.09866367e-01;
  rot_mat[3][3] = 5.55699395e-01;
  rot_mat[3][4] = 5.73062990e-01;
  rot_mat[3][5] = 1.74778224e-01;
  rot_mat[3][6] = 8.77170506e-03;
  rot_mat[4][0] = -3.24609748e-01;
  rot_mat[4][1] = 2.89664179e-01;
  rot_mat[4][2] = -7.50613752e-01;
  rot_mat[4][3] = -1.63060005e-01;
  rot_mat[4][4] = 1.70803377e-01;
  rot_mat[4][5] = 1.63784167e-01;
  rot_mat[4][6] = 4.05850414e-01;
  rot_mat[5][0] = 1.32570771e-03;
  rot_mat[5][1] = 3.80299846e-01;
  rot_mat[5][2] = -5.08439810e-02;
  rot_mat[5][3] = 7.59141791e-01;
  rot_mat[5][4] = -4.77844928e-01;
  rot_mat[5][5] = 2.12149087e-01;
  rot_mat[5][6] = 5.60882349e-02;
  rot_mat[6][0] = 1.62781208e-01;
  rot_mat[6][1] = -5.75073150e-01;
  rot_mat[6][2] = -3.60485665e-01;
  rot_mat[6][3] = -1.70070331e-02;
  rot_mat[6][4] = -5.72251258e-01;
  rot_mat[6][5] = -3.50549638e-01;
  rot_mat[6][6] = 2.49394158e-01;
  spo->storeParamsBeforeRotation(); // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat, false);

  // Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Compute the expected value, gradient, and laplacian by hand using the transformation above
  // Check value
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        psiM_rot_manual[i][j] += psiM_bare[i][k] * rot_mat[k][j];
    }
  checkMatrix(psiM_rot_manual, psiM_rot);

  // Check grad
  SPOSet::GradMatrix dpsiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      dpsiM_rot_manual[i][j][0] = 0.;
      dpsiM_rot_manual[i][j][1] = 0.;
      dpsiM_rot_manual[i][j][2] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        for (int l = 0; l < 3; l++)
          dpsiM_rot_manual[i][j][l] += dpsiM_bare[i][k][l] * rot_mat[k][j];
    }

  // No checkMatrix for tensors? Gotta do it by hand...
  double res = 0.;
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
      for (int k = 0; k < 3; k++)
        res += std::sqrt(std::norm(dpsiM_rot_manual[i][j][k] - dpsiM_rot[i][j][k]));

  CHECK(res == Approx(0.).epsilon(2e-4)); // increase threshold to accomodate mixed precision.

  // Check laplacian
  SPOSet::ValueMatrix d2psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      d2psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        d2psiM_rot_manual[i][j] += d2psiM_bare[i][k] * rot_mat[k][j];
    }

  checkMatrix(d2psiM_rot_manual, d2psiM_rot);

} // TEST_CASE


TEST_CASE("Spline applyRotation two rotations", "[wavefunction]")
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
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" gpu="no" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Check before rotation is applied. From 'test_einset.cpp'
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

  /* Apply a rotation, U = exp(-K) as shown below, to the splines.
     Because K = -K.T, U is unitary.
     K=
      0.00000000e+00  -2.33449314e-01   5.21735139e-01   1.67744276e-01  -1.68565994e-01   1.26632312e-02   2.29272040e-01
      2.33449314e-01   0.00000000e+00  -3.81982621e-01  -6.76995896e-01   7.15727948e-02   9.10926961e-01  -6.94864205e-01
     -5.21735139e-01   3.81982621e-01   0.00000000e+00   2.37433139e-01  -7.09878105e-01  -6.33841289e-01  -7.49009582e-01
     -1.67744276e-01   6.76995896e-01  -2.37433139e-01   0.00000000e+00  -6.72002172e-01   4.27152921e-01  -4.82983743e-03
      1.68565994e-01  -7.15727948e-02   7.09878105e-01   6.72002172e-01   0.00000000e+00  -7.93860012e-01  -7.76624731e-01
     -1.26632312e-02  -9.10926961e-01   6.33841289e-01  -4.27152921e-01   7.93860012e-01   0.00000000e+00  -2.74931052e-01
     -2.29272040e-01   6.94864205e-01   7.49009582e-01   4.82983743e-03   7.76624731e-01   2.74931052e-01   0.00000000e+00
     U=
      8.06061880e-01   3.54921598e-01  -3.40706426e-01  -5.90163619e-02   1.11650454e-01  -1.99768450e-01  -2.28818375e-01
      1.58069821e-02   2.10363421e-01   2.74922448e-01   2.12581764e-01   2.64602356e-01  -6.34971914e-01   6.01265560e-01
      4.43696646e-01  -6.06912539e-02   2.61413193e-01  -1.98368802e-01  -6.43234645e-02   5.75880430e-01   5.96646495e-01
      1.45865363e-01  -5.16577220e-01  -2.09866367e-01   5.55699395e-01   5.73062990e-01   1.74778224e-01   8.77170506e-03
     -3.24609748e-01   2.89664179e-01  -7.50613752e-01  -1.63060005e-01   1.70803377e-01   1.63784167e-01   4.05850414e-01
      1.32570771e-03   3.80299846e-01  -5.08439810e-02   7.59141791e-01  -4.77844928e-01   2.12149087e-01   5.60882349e-02
      1.62781208e-01  -5.75073150e-01  -3.60485665e-01  -1.70070331e-02  -5.72251258e-01  -3.50549638e-01   2.49394158e-01
     UU.T=
      1.00000000e+00  -2.77555756e-17   5.55111512e-17  -8.67361738e-18  -1.11022302e-16  -6.07153217e-17   6.93889390e-17
     -2.77555756e-17   1.00000000e+00   5.55111512e-17  -1.84748050e-16   5.55111512e-17  -6.93889390e-18   8.32667268e-17
      5.55111512e-17   5.55111512e-17   1.00000000e+00   1.45716772e-16  -5.55111512e-17   1.38777878e-16   5.55111512e-17
     -8.67361738e-18  -1.84748050e-16   1.45716772e-16   1.00000000e+00  -1.95156391e-17   1.80411242e-16   8.02309608e-17
     -1.11022302e-16   5.55111512e-17  -5.55111512e-17  -1.95156391e-17   1.00000000e+00   7.28583860e-17   1.11022302e-16
     -6.07153217e-17  -6.93889390e-18   1.38777878e-16   1.80411242e-16   7.28583860e-17   1.00000000e+00  -1.12757026e-16
      6.93889390e-17   8.32667268e-17   5.55111512e-17   8.02309608e-17   1.11022302e-16  -1.12757026e-16   1.00000000e+00

      NB: There's nothing special about this rotation. I purposefully made something 'ugly' to make a worst case scenario...
    */
  SPOSet::ValueMatrix rot_mat1(orbitalsetsize, orbitalsetsize);
  rot_mat1[0][0] = 8.06061880e-01;
  rot_mat1[0][1] = 3.54921598e-01;
  rot_mat1[0][2] = -3.40706426e-01;
  rot_mat1[0][3] = -5.90163619e-02;
  rot_mat1[0][4] = 1.11650454e-01;
  rot_mat1[0][5] = -1.99768450e-01;
  rot_mat1[0][6] = -2.28818375e-01;
  rot_mat1[1][0] = 1.58069821e-02;
  rot_mat1[1][1] = 2.10363421e-01;
  rot_mat1[1][2] = 2.74922448e-01;
  rot_mat1[1][3] = 2.12581764e-01;
  rot_mat1[1][4] = 2.64602356e-01;
  rot_mat1[1][5] = -6.34971914e-01;
  rot_mat1[1][6] = 6.01265560e-01;
  rot_mat1[2][0] = 4.43696646e-01;
  rot_mat1[2][1] = -6.06912539e-02;
  rot_mat1[2][2] = 2.61413193e-01;
  rot_mat1[2][3] = -1.98368802e-01;
  rot_mat1[2][4] = -6.43234645e-02;
  rot_mat1[2][5] = 5.75880430e-01;
  rot_mat1[2][6] = 5.96646495e-01;
  rot_mat1[3][0] = 1.45865363e-01;
  rot_mat1[3][1] = -5.16577220e-01;
  rot_mat1[3][2] = -2.09866367e-01;
  rot_mat1[3][3] = 5.55699395e-01;
  rot_mat1[3][4] = 5.73062990e-01;
  rot_mat1[3][5] = 1.74778224e-01;
  rot_mat1[3][6] = 8.77170506e-03;
  rot_mat1[4][0] = -3.24609748e-01;
  rot_mat1[4][1] = 2.89664179e-01;
  rot_mat1[4][2] = -7.50613752e-01;
  rot_mat1[4][3] = -1.63060005e-01;
  rot_mat1[4][4] = 1.70803377e-01;
  rot_mat1[4][5] = 1.63784167e-01;
  rot_mat1[4][6] = 4.05850414e-01;
  rot_mat1[5][0] = 1.32570771e-03;
  rot_mat1[5][1] = 3.80299846e-01;
  rot_mat1[5][2] = -5.08439810e-02;
  rot_mat1[5][3] = 7.59141791e-01;
  rot_mat1[5][4] = -4.77844928e-01;
  rot_mat1[5][5] = 2.12149087e-01;
  rot_mat1[5][6] = 5.60882349e-02;
  rot_mat1[6][0] = 1.62781208e-01;
  rot_mat1[6][1] = -5.75073150e-01;
  rot_mat1[6][2] = -3.60485665e-01;
  rot_mat1[6][3] = -1.70070331e-02;
  rot_mat1[6][4] = -5.72251258e-01;
  rot_mat1[6][5] = -3.50549638e-01;
  rot_mat1[6][6] = 2.49394158e-01;
  spo->storeParamsBeforeRotation();    // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat1, false); // Apply rotation 1


  /* Apply another rotation, U = exp(-K) as shown below, to the splines.
     Because K = -K.T, U is unitary.
     K=
      0.00000000e+00   3.08898211e-01   9.70811450e-01   9.96582440e-01   2.59290113e-01   9.08544511e-01   7.47970513e-01
     -3.08898211e-01   0.00000000e+00   4.90958419e-01   7.98394113e-01   9.02926177e-01   3.24156332e-01   1.83039207e-01
     -9.70811450e-01  -4.90958419e-01   0.00000000e+00   7.73447329e-01   2.71156433e-01   2.18009012e-01   1.75304629e-01
     -9.96582440e-01  -7.98394113e-01  -7.73447329e-01   0.00000000e+00   9.27679516e-01   8.09853231e-01   7.71485500e-01
     -2.59290113e-01  -9.02926177e-01  -2.71156433e-01  -9.27679516e-01   0.00000000e+00   6.19032776e-01   6.05979744e-02
     -9.08544511e-01  -3.24156332e-01  -2.18009012e-01  -8.09853231e-01  -6.19032776e-01   0.00000000e+00   7.91708106e-01
     -7.47970513e-01  -1.83039207e-01  -1.75304629e-01  -7.71485500e-01  -6.05979744e-02  -7.91708106e-01   0.00000000e+00

     U=
     -2.98194829e-02  -6.14346084e-01  -6.43923495e-01  -1.18552217e-01   3.77244445e-01  -2.06047353e-01   9.07122365e-02
     -3.48818370e-01   4.38506143e-01  -5.57140467e-01  -3.70956624e-01  -2.23467853e-01   3.80270838e-01   2.08518583e-01
      1.87644050e-01  -1.34182635e-01   3.74618322e-01  -8.74132032e-01   1.57504581e-01   4.78555353e-02   1.23455215e-01
     -7.68247448e-02  -2.48006571e-01  -7.30866440e-02  -2.06817660e-01  -7.99320897e-01  -4.53458397e-01  -1.99842648e-01
      2.59360916e-02   5.64417889e-01  -1.05874452e-01  -1.09701962e-01   2.76429977e-01  -7.59899812e-01   6.04531910e-02
      3.50608068e-01   1.58922313e-01  -2.49914906e-01  -1.49783413e-01   1.03865829e-01   1.56180926e-01  -8.55420691e-01
      8.44230723e-01   8.33975900e-02  -2.35816939e-01   8.35859456e-02  -2.38381031e-01   5.64243083e-02   3.97132056e-01

     UU.T=
      1.00000000e+00   1.97376302e-16   2.32170398e-16   8.48209595e-17   1.52095225e-16   1.20549468e-16  -6.20012223e-17
      1.97376302e-16   1.00000000e+00   2.49491048e-16  -1.32824448e-16  -2.95131454e-17   1.89706297e-17  -9.33948863e-17
      2.32170398e-16   2.49491048e-16   1.00000000e+00  -1.53719614e-16  -1.07440039e-16   1.11845140e-17   1.25992152e-16
      8.48209595e-17  -1.32824448e-16  -1.53719614e-16   1.00000000e+00   3.98945291e-17   8.97184517e-17  -1.23231760e-16
      1.52095225e-16  -2.95131454e-17  -1.07440039e-16   3.98945291e-17   1.00000000e+00   2.40723889e-17   3.26140430e-17
      1.20549468e-16   1.89706297e-17   1.11845140e-17   8.97184517e-17   2.40723889e-17   1.00000000e+00   2.75131978e-17
     -6.20012223e-17  -9.33948863e-17   1.25992152e-16  -1.23231760e-16   3.26140430e-17   2.75131978e-17   1.00000000e+00

      NB: There's nothing special about this rotation. I purposefully made something 'ugly' to make a worst case scenario...
    */
  SPOSet::ValueMatrix rot_mat2(orbitalsetsize, orbitalsetsize);
  rot_mat2[0][0] = -2.98194829e-02;
  rot_mat2[0][1] = -6.14346084e-01;
  rot_mat2[0][2] = -6.43923495e-01;
  rot_mat2[0][3] = -1.18552217e-01;
  rot_mat2[0][4] = 3.77244445e-01;
  rot_mat2[0][5] = -2.06047353e-01;
  rot_mat2[0][6] = 9.07122365e-02;
  rot_mat2[1][0] = -3.48818370e-01;
  rot_mat2[1][1] = 4.38506143e-01;
  rot_mat2[1][2] = -5.57140467e-01;
  rot_mat2[1][3] = -3.70956624e-01;
  rot_mat2[1][4] = -2.23467853e-01;
  rot_mat2[1][5] = 3.80270838e-01;
  rot_mat2[1][6] = 2.08518583e-01;
  rot_mat2[2][0] = 1.87644050e-01;
  rot_mat2[2][1] = -1.34182635e-01;
  rot_mat2[2][2] = 3.74618322e-01;
  rot_mat2[2][3] = -8.74132032e-01;
  rot_mat2[2][4] = 1.57504581e-01;
  rot_mat2[2][5] = 4.78555353e-02;
  rot_mat2[2][6] = 1.23455215e-01;
  rot_mat2[3][0] = -7.68247448e-02;
  rot_mat2[3][1] = -2.48006571e-01;
  rot_mat2[3][2] = -7.30866440e-02;
  rot_mat2[3][3] = -2.06817660e-01;
  rot_mat2[3][4] = -7.99320897e-01;
  rot_mat2[3][5] = -4.53458397e-01;
  rot_mat2[3][6] = -1.99842648e-01;
  rot_mat2[4][0] = 2.59360916e-02;
  rot_mat2[4][1] = 5.64417889e-01;
  rot_mat2[4][2] = -1.05874452e-01;
  rot_mat2[4][3] = -1.09701962e-01;
  rot_mat2[4][4] = 2.76429977e-01;
  rot_mat2[4][5] = -7.59899812e-01;
  rot_mat2[4][6] = 6.04531910e-02;
  rot_mat2[5][0] = 3.50608068e-01;
  rot_mat2[5][1] = 1.58922313e-01;
  rot_mat2[5][2] = -2.49914906e-01;
  rot_mat2[5][3] = -1.49783413e-01;
  rot_mat2[5][4] = 1.03865829e-01;
  rot_mat2[5][5] = 1.56180926e-01;
  rot_mat2[5][6] = -8.55420691e-01;
  rot_mat2[6][0] = 8.44230723e-01;
  rot_mat2[6][1] = 8.33975900e-02;
  rot_mat2[6][2] = -2.35816939e-01;
  rot_mat2[6][3] = 8.35859456e-02;
  rot_mat2[6][4] = -2.38381031e-01;
  rot_mat2[6][5] = 5.64243083e-02;
  rot_mat2[6][6] = 3.97132056e-01;
  spo->storeParamsBeforeRotation();    // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat2, false); // Apply rotation2

  // Total rotation is product of rot_mat1 and rot_mat2
  SPOSet::ValueMatrix rot_mat_tot(orbitalsetsize, orbitalsetsize);
  for (int i = 0; i < orbitalsetsize; i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      rot_mat_tot[i][j] = 0;
      for (int k = 0; k < orbitalsetsize; k++)
        rot_mat_tot[i][j] += rot_mat1[i][k] * rot_mat2[k][j];
    }

  // Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Compute the expected value, gradient, and laplacian by hand using the transformation above
  // Check value
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        psiM_rot_manual[i][j] += psiM_bare[i][k] * rot_mat_tot[k][j];
    }
  checkMatrix(psiM_rot_manual, psiM_rot);

  // Check grad
  SPOSet::GradMatrix dpsiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      dpsiM_rot_manual[i][j][0] = 0.;
      dpsiM_rot_manual[i][j][1] = 0.;
      dpsiM_rot_manual[i][j][2] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        for (int l = 0; l < 3; l++)
          dpsiM_rot_manual[i][j][l] += dpsiM_bare[i][k][l] * rot_mat_tot[k][j];
    }

  // No checkMatrix for tensors? Gotta do it by hand...
  double res = 0.;
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
      for (int k = 0; k < 3; k++)
        res += std::sqrt(std::norm(dpsiM_rot_manual[i][j][k] - dpsiM_rot[i][j][k]));

  CHECK(res == Approx(0.).epsilon(2e-4)); // increase threshold to accomodate mixed precision.

  // Check laplacian
  SPOSet::ValueMatrix d2psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      d2psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        d2psiM_rot_manual[i][j] += d2psiM_bare[i][k] * rot_mat_tot[k][j];
    }
  checkMatrix(d2psiM_rot_manual, d2psiM_rot);

} // TEST_CASE


#ifdef QMC_COMPLEX
TEST_CASE("Spline applyRotation complex rotation", "[wavefunction]")
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
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" gpu="no" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Check before rotation is applied. From 'test_einset.cpp'
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

  /* Apply a rotation, U = exp(-K) as shown below, to the splines.
     Because K = -K.T, U is unitary.
     K= (0+i)*eye(orbitalsetsize)
     U= (5.40302306e-01 -8.41470985e-01i) * eye(orbitalsetsize)
     U * conj(U.T) = eye(7)

     NB: There's nothing special about this rotation. I purposefully made something 'ugly' to make a worst case scenario...
    */
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  const SPOSet::ValueType z(5.40302306e-01, -8.41470985e-01);
  rot_mat = 0;
  for (int i = 0; i < orbitalsetsize; i++)
    rot_mat[i][i] = z;

  spo->storeParamsBeforeRotation(); // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat, false);

  // Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Compute the expected value, gradient, and laplacian by hand using the transformation above
  // Check value
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        psiM_rot_manual[i][j] += psiM_bare[i][k] * rot_mat[k][j];
    }
  checkMatrix(psiM_rot_manual, psiM_rot);

  // Check grad
  SPOSet::GradMatrix dpsiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      dpsiM_rot_manual[i][j][0] = 0.;
      dpsiM_rot_manual[i][j][1] = 0.;
      dpsiM_rot_manual[i][j][2] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        for (int l = 0; l < 3; l++)
          dpsiM_rot_manual[i][j][l] += dpsiM_bare[i][k][l] * rot_mat[k][j];
    }

  // No checkMatrix for tensors? Gotta do it by hand...
  double res = 0.;
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
      for (int k = 0; k < 3; k++)
        res += std::sqrt(std::norm(dpsiM_rot_manual[i][j][k] - dpsiM_rot[i][j][k]));

  CHECK(res == Approx(0.).epsilon(2e-4)); // increase threshold to accomodate mixed precision.

  // Check laplacian
  SPOSet::ValueMatrix d2psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (int i = 0; i < elec_.R.size(); i++)
    for (int j = 0; j < orbitalsetsize; j++)
    {
      d2psiM_rot_manual[i][j] = 0.;
      for (int k = 0; k < orbitalsetsize; k++)
        d2psiM_rot_manual[i][j] += d2psiM_bare[i][k] * rot_mat[k][j];
    }
  checkMatrix(d2psiM_rot_manual, d2psiM_rot);
} // TEST_CASE
#endif
} // namespace qmcplusplus
