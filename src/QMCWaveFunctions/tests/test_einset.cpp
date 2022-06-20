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
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Einspline SPO from HDF diamond_1x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // monoO
  /*
  lattice.R(0,0) = 5.10509515;
  lattice.R(0,1) = -3.23993545;
  lattice.R(0,2) = 0.0;
  lattice.R(1,0) = 5.10509515;
  lattice.R(1,1) = 3.23993545;
  lattice.R(1,2) = 0.0;
  lattice.R(2,0) = -6.49690625;
  lattice.R(2,1) = 0.0;
  lattice.R(2,2) = 7.08268015;
  */

  // diamondC_1x1x1
  lattice.R(0, 0) = 3.37316115;
  lattice.R(0, 1) = 3.37316115;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.37316115;
  lattice.R(1, 2) = 3.37316115;
  lattice.R(2, 0) = 3.37316115;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.37316115;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0][0] = 0.0;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"8\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // value
  REQUIRE(std::real(psiM[1][0]) == Approx(-0.8886948824));
  REQUIRE(std::real(psiM[1][1]) == Approx(1.4194120169));
  // grad
  REQUIRE(std::real(dpsiM[1][0][0]) == Approx(-0.0000183403));
  REQUIRE(std::real(dpsiM[1][0][1]) == Approx(0.1655139178));
  REQUIRE(std::real(dpsiM[1][0][2]) == Approx(-0.0000193077));
  REQUIRE(std::real(dpsiM[1][1][0]) == Approx(-1.3131694794));
  REQUIRE(std::real(dpsiM[1][1][1]) == Approx(-1.1174004078));
  REQUIRE(std::real(dpsiM[1][1][2]) == Approx(-0.8462534547));
  // lapl
  REQUIRE(std::real(d2psiM[1][0]) == Approx(1.3313053846));
  REQUIRE(std::real(d2psiM[1][1]) == Approx(-4.712583065));

  // for vgh
  SPOSet::ValueVector psiV(psiM[1], spo->getOrbitalSetSize());
  SPOSet::GradVector dpsiV(dpsiM[1], spo->getOrbitalSetSize());
  SPOSet::HessVector ddpsiV(spo->getOrbitalSetSize());
  spo->evaluateVGH(elec_, 1, psiV, dpsiV, ddpsiV);

  // Catch default is 100*(float epsilson)
  double eps = 2000 * std::numeric_limits<float>::epsilon();
  //hess
  REQUIRE(std::real(ddpsiV[1](0, 0)) == Approx(-2.3160984034));
  REQUIRE(std::real(ddpsiV[1](0, 1)) == Approx(1.8089479397));
  REQUIRE(std::real(ddpsiV[1](0, 2)) == Approx(0.5608575749));
  REQUIRE(std::real(ddpsiV[1](1, 0)) == Approx(1.8089479397));
  REQUIRE(std::real(ddpsiV[1](1, 1)) == Approx(-0.07996207476).epsilon(eps));
  REQUIRE(std::real(ddpsiV[1](1, 2)) == Approx(0.5237969314));
  REQUIRE(std::real(ddpsiV[1](2, 0)) == Approx(0.5608575749));
  REQUIRE(std::real(ddpsiV[1](2, 1)) == Approx(0.5237969314));
  REQUIRE(std::real(ddpsiV[1](2, 2)) == Approx(-2.316497764));

  SPOSet::HessMatrix hesspsiV(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GGGMatrix d3psiV(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, hesspsiV, d3psiV);

  //The reference values for grad_grad_grad_psi.
  /*
  d3psiV(1,0)[0][0]=(0.046337127685546875,-0.046337127685546875)
  d3psiV(1,0)[0][1]=(1.1755813360214233,-1.1755813360214233)
  d3psiV(1,0)[0][2]=(0.066015571355819702,-0.066015541553497314)
  d3psiV(1,0)[0][4]=(0.041470438241958618,-0.041470438241958618)
  d3psiV(1,0)[0][5]=(-0.51674127578735352,0.51674127578735352)
  d3psiV(1,0)[0][8]=(0.065953642129898071,-0.065953642129898071)
  d3psiV(1,0)[1][4]=(-4.8771157264709473,4.8771157264709473)
  d3psiV(1,0)[1][5]=(0.041532635688781738,-0.041532605886459351)
  d3psiV(1,0)[1][8]=(1.1755810976028442,-1.1755810976028442)
  d3psiV(1,0)[2][8]=(0.046399354934692383,-0.046399354934692383)

  d3psiV(1,1)[0][0]=(6.7155771255493164,-7.5906991958618164)
  d3psiV(1,1)[0][1]=(5.545051097869873,-5.0280308723449707)
  d3psiV(1,1)[0][2]=(0.98297119140625,-0.50021600723266602)
  d3psiV(1,1)[0][4]=(-3.1704092025756836,3.8900821208953857)
  d3psiV(1,1)[0][5]=(-1.9537661075592041,1.7758266925811768)
  d3psiV(1,1)[0][8]=(1.9305641651153564,-2.1480715274810791)
  d3psiV(1,1)[1][4]=(3.605137825012207,-3.2767453193664551)
  d3psiV(1,1)[1][5]=(-0.73825764656066895,-0.33745908737182617)
  d3psiV(1,1)[1][8]=(5.5741839408874512,-5.0784988403320312)
  d3psiV(1,1)[2][8]=(3.131234884262085,-1.3596141338348389)
*/

#if 0 //Enable when finite precision issue on Rhea is found.
  REQUIRE(std::real(d3psiV(1,0)[0][0] ) == Approx(0.0463371276));
  REQUIRE(std::real(d3psiV(1,0)[0][1] ) == Approx(1.1755813360));
  REQUIRE(std::real(d3psiV(1,0)[0][2] ) == Approx(0.0660155713));
  REQUIRE(std::real(d3psiV(1,0)[0][4] ) == Approx(0.0414704382));
  REQUIRE(std::real(d3psiV(1,0)[0][5] ) == Approx(-0.5167412758));
  REQUIRE(std::real(d3psiV(1,0)[0][8] ) == Approx(0.0659536421));
  REQUIRE(std::real(d3psiV(1,0)[1][4] ) == Approx(-4.8771157264));
  REQUIRE(std::real(d3psiV(1,0)[1][5] ) == Approx(0.0415326356));
  REQUIRE(std::real(d3psiV(1,0)[1][8] ) == Approx(1.1755810976));
  REQUIRE(std::real(d3psiV(1,0)[2][8] ) == Approx(0.0463993549));

  REQUIRE(std::real(d3psiV(1,1)[0][0] ) == Approx(6.7155771255));
  REQUIRE(std::real(d3psiV(1,1)[0][1] ) == Approx(5.5450510978));
  REQUIRE(std::real(d3psiV(1,1)[0][2] ) == Approx(0.9829711914));
  REQUIRE(std::real(d3psiV(1,1)[0][4] ) == Approx(-3.1704092025));
  REQUIRE(std::real(d3psiV(1,1)[0][5] ) == Approx(-1.9537661075));
  REQUIRE(std::real(d3psiV(1,1)[0][8] ) == Approx(1.9305641651));
  REQUIRE(std::real(d3psiV(1,1)[1][4] ) == Approx(3.6051378250));
  REQUIRE(std::real(d3psiV(1,1)[1][5] ) == Approx(-0.7382576465));
  REQUIRE(std::real(d3psiV(1,1)[1][8] ) == Approx(5.5741839408));
  REQUIRE(std::real(d3psiV(1,1)[2][8] ) == Approx(3.1312348842));
#endif

#endif

  // SplineR2R only for the moment, so skip if QMC_COMPLEX is set
#if !defined(QMC_CUDA) && !defined(QMC_COMPLEX)
  /* 
     Here we test the low-level spline implementation of applyRotation().

     Manually encode the unitary transformation. Ugly, but it works.
     @TODO: Use the total rotation machinery when it's implemented

     NB: This is truncated to 5 sig-figs, so there is some slop here as
         compared to what is done in the splines via apply_rotation().
	 So below we reduce the threshold for comparison. This can
	 probably be ditched once we have a way to grab the actual 
	 rotation matrix...
  */
  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 8);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // value
  CHECK(psiM_bare[1][0] == Approx(-0.8886948824));
  CHECK(psiM_bare[1][1] == Approx(1.4194120169));
  // grad
  CHECK(dpsiM_bare[1][0][0] == Approx(-0.0000183403));
  CHECK(dpsiM_bare[1][0][1] == Approx(0.1655139178));
  CHECK(dpsiM_bare[1][0][2] == Approx(-0.0000193077));
  CHECK(dpsiM_bare[1][1][0] == Approx(-1.3131694794));
  CHECK(dpsiM_bare[1][1][1] == Approx(-1.1174004078));
  CHECK(dpsiM_bare[1][1][2] == Approx(-0.8462534547));
  // lapl
  CHECK(d2psiM_bare[1][0] == Approx(1.3313053846));
  CHECK(d2psiM_bare[1][1] == Approx(-4.712583065));

  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat[0][0] = 0.99726;
  rot_mat[0][1] = -0.00722;
  rot_mat[0][2] = 0.00014;
  rot_mat[0][3] = -0.00982;
  rot_mat[0][4] = -0.01979;
  rot_mat[0][5] = -0.02976;
  rot_mat[0][6] = -0.03972;
  rot_mat[0][7] = -0.04969;
  rot_mat[1][0] = -0.00722;
  rot_mat[1][1] = 0.97754;
  rot_mat[1][2] = -0.05955;
  rot_mat[1][3] = -0.06945;
  rot_mat[1][4] = -0.07935;
  rot_mat[1][5] = -0.08925;
  rot_mat[1][6] = -0.09915;
  rot_mat[1][7] = -0.10905;
  rot_mat[2][0] = -0.00014;
  rot_mat[2][1] = 0.05955;
  rot_mat[2][2] = 0.99821;
  rot_mat[2][3] = -0.00209;
  rot_mat[2][4] = -0.00239;
  rot_mat[2][5] = -0.00269;
  rot_mat[2][6] = -0.00299;
  rot_mat[2][7] = -0.00329;
  rot_mat[3][0] = 0.00982;
  rot_mat[3][1] = 0.06945;
  rot_mat[3][2] = -0.00209;
  rot_mat[3][3] = 0.99751;
  rot_mat[3][4] = -0.00289;
  rot_mat[3][5] = -0.00329;
  rot_mat[3][6] = -0.00368;
  rot_mat[3][7] = -0.00408;
  rot_mat[4][0] = 0.01979;
  rot_mat[4][1] = 0.07935;
  rot_mat[4][2] = -0.00239;
  rot_mat[4][3] = -0.00289;
  rot_mat[4][4] = 0.99661;
  rot_mat[4][5] = -0.00388;
  rot_mat[4][6] = -0.00438;
  rot_mat[4][7] = -0.00488;
  rot_mat[5][0] = 0.02976;
  rot_mat[5][1] = 0.08925;
  rot_mat[5][2] = -0.00269;
  rot_mat[5][3] = -0.00329;
  rot_mat[5][4] = -0.00388;
  rot_mat[5][5] = 0.99552;
  rot_mat[5][6] = -0.00508;
  rot_mat[5][7] = -0.00568;
  rot_mat[6][0] = 0.03972;
  rot_mat[6][1] = 0.09915;
  rot_mat[6][2] = -0.00299;
  rot_mat[6][3] = -0.00368;
  rot_mat[6][4] = -0.00438;
  rot_mat[6][5] = -0.00508;
  rot_mat[6][6] = 0.99422;
  rot_mat[6][7] = -0.00647;
  rot_mat[7][0] = 0.04969;
  rot_mat[7][1] = 0.10905;
  rot_mat[7][2] = -0.00329;
  rot_mat[7][3] = -0.00408;
  rot_mat[7][4] = -0.00488;
  rot_mat[7][5] = -0.00568;
  rot_mat[7][6] = -0.00647;
  rot_mat[7][7] = 0.99273;

  // Apply the rotation
  spo->applyRotation(rot_mat, false);

  // Get data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Now compute the expected values by hand using the transformation above
  double val1 = 0.;
  double val2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    val1 += psiM_bare[0][i] * rot_mat[i][0];
    val2 += psiM_bare[1][i] * rot_mat[i][0];
  }

  // value
  CHECK(psiM_rot[0][0] == Approx(val1));
  CHECK(psiM_rot[1][0] == Approx(val2));

  std::vector<double> grad1(3);
  std::vector<double> grad2(3);
  for (auto j = 0; j < grad1.size(); j++)
  {
    for (auto i = 0; i < rot_mat.size1(); i++)
    {
      grad1[j] += dpsiM_bare[0][i][j] * rot_mat[i][0];
      grad2[j] += dpsiM_bare[1][i][j] * rot_mat[i][0];
    }
  }

  // grad
  CHECK(dpsiM_rot[0][0][0] == Approx(grad1[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][1] == Approx(grad1[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][2] == Approx(grad1[2]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][0] == Approx(grad2[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][1] == Approx(grad2[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][2] == Approx(grad2[2]).epsilon(0.0001));

  double lap1 = 0.;
  double lap2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    lap1 += d2psiM_bare[0][i] * rot_mat[i][0];
    lap2 += d2psiM_bare[1][i] * rot_mat[i][0];
  }

  // Lapl
  CHECK(d2psiM_rot[0][0] == Approx(lap1).epsilon(0.0001));
  CHECK(d2psiM_rot[1][0] == Approx(lap2).epsilon(0.0001));
#endif

#if 0
  // Dump values of the orbitals
  int orbSize= spo->getOrbitalSetSize();
  int basisSize= spo->getBasisSetSize();
  printf("orb size = %d basis set size = %d\n",orbSize, basisSize);

  FILE *fspo = fopen("spo.dat", "w");
  for (int ix = 0; ix < 30; ix++) {
    for (int iy = 0; iy < 30; iy++) {
      for (int iz = 0; iz < 30; iz++) {
        double x = 0.1*ix - 1.5;
        double y = 0.1*iy - 1.5;
        double z = 0.1*iz - 1.5;
        elec_.R[0][0] = x;
        elec_.R[0][1] = y;
        elec_.R[0][2] = z;
        elec_.update();
        SPOSet::ValueVector orbs(orbSize);
        spo->evaluate(elec_, 0, orbs);
        fprintf(fspo, "%g %g %g",x,y,z);
        for (int j = 0; j < orbSize; j++) {
          fprintf(fspo, " %g ",orbs[j]);
        }
        fprintf(fspo, "\n");
      }
    }
  }
  fclose(fspo);

#endif
}

TEST_CASE("Einspline SPO from HDF diamond_2x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // diamondC_2x1x1
  lattice.R(0, 0) = 6.7463223;
  lattice.R(0, 1) = 6.7463223;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.37316115;
  lattice.R(1, 2) = 3.37316115;
  lattice.R(2, 0) = 3.37316115;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.37316115;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({4});
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;
  ions_.R[2][0] = 3.37316115;
  ions_.R[2][1] = 3.37316115;
  ions_.R[2][2] = 0.0;
  ions_.R[3][0] = 5.05974173;
  ions_.R[3][1] = 5.05974173;
  ions_.R[3][2] = 1.68658058;


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0][0] = 0.0;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_2x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_2x1x1.pwscf.h5\" tilematrix=\"2 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  // for vgl
  SPOSet::ValueMatrix psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  REQUIRE(std::real(psiM[1][0]) == Approx(0.9008999467));
  REQUIRE(std::real(psiM[1][1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::real(dpsiM[1][0][0]) == Approx(0.0025820041));
  REQUIRE(std::real(dpsiM[1][0][1]) == Approx(-0.1880052537));
  REQUIRE(std::real(dpsiM[1][0][2]) == Approx(-0.0025404284));
  REQUIRE(std::real(dpsiM[1][1][0]) == Approx(0.1069662273));
  REQUIRE(std::real(dpsiM[1][1][1]) == Approx(-0.4364597797));
  REQUIRE(std::real(dpsiM[1][1][2]) == Approx(-0.106951952));
  // lapl
  REQUIRE(std::real(d2psiM[1][0]) == Approx(-1.3757134676));
  REQUIRE(std::real(d2psiM[1][1]) == Approx(-2.4803137779));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  REQUIRE(std::imag(psiM[1][0]) == Approx(0.9008999467));
  REQUIRE(std::imag(psiM[1][1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::imag(dpsiM[1][0][0]) == Approx(0.0025820041));
  REQUIRE(std::imag(dpsiM[1][0][1]) == Approx(-0.1880052537));
  REQUIRE(std::imag(dpsiM[1][0][2]) == Approx(-0.0025404284));
  REQUIRE(std::imag(dpsiM[1][1][0]) == Approx(0.1069453433));
  REQUIRE(std::imag(dpsiM[1][1][1]) == Approx(-0.43649593));
  REQUIRE(std::imag(dpsiM[1][1][2]) == Approx(-0.1069145575));
  // lapl
  REQUIRE(std::imag(d2psiM[1][0]) == Approx(-1.3757134676));
  REQUIRE(std::imag(d2psiM[1][1]) == Approx(-2.4919104576));
#endif

  // test batched interfaces

  ParticleSet elec_2(elec_);
  // interchange positions
  elec_2.R[0] = elec_.R[1];
  elec_2.R[1] = elec_.R[0];
  RefVectorWithLeader<ParticleSet> p_list(elec_);
  p_list.push_back(elec_);
  p_list.push_back(elec_2);

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  RefVectorWithLeader<SPOSet> spo_list(*spo);
  spo_list.push_back(*spo);
  spo_list.push_back(*spo_2);

  SPOSet::ValueVector psi(spo->getOrbitalSetSize());
  SPOSet::GradVector dpsi(spo->getOrbitalSetSize());
  SPOSet::ValueVector d2psi(spo->getOrbitalSetSize());
  SPOSet::ValueVector psi_2(spo->getOrbitalSetSize());
  SPOSet::GradVector dpsi_2(spo->getOrbitalSetSize());
  SPOSet::ValueVector d2psi_2(spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueVector> psi_v_list;
  RefVector<SPOSet::GradVector> dpsi_v_list;
  RefVector<SPOSet::ValueVector> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);

  spo->mw_evaluateVGL(spo_list, p_list, 0, psi_v_list, dpsi_v_list, d2psi_v_list);
#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  REQUIRE(std::real(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::real(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::real(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  REQUIRE(std::real(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  REQUIRE(std::real(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][0]) == Approx(0.1069662273));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][1]) == Approx(-0.4364597797));
  REQUIRE(std::real(dpsi_v_list[1].get()[1][2]) == Approx(-0.106951952));
  // lapl
  REQUIRE(std::real(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  REQUIRE(std::real(d2psi_v_list[1].get()[1]) == Approx(-2.4803137779));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  REQUIRE(std::imag(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::imag(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  REQUIRE(std::imag(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][0]) == Approx(0.1069453433));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][1]) == Approx(-0.43649593));
  REQUIRE(std::imag(dpsi_v_list[1].get()[1][2]) == Approx(-0.1069145575));
  // lapl
  REQUIRE(std::imag(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  REQUIRE(std::imag(d2psi_v_list[1].get()[1]) == Approx(-2.4919104576));
#endif
}

TEST_CASE("EinsplineSetBuilder CheckLattice", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  lattice.R       = 0.0;
  lattice.R(0, 0) = 1.0;
  lattice.R(1, 1) = 1.0;
  lattice.R(2, 2) = 1.0;

  const SimulationCell simulation_cell(lattice);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto& elec(*elec_ptr);

  elec.setName("elec");
  std::vector<int> agroup(2);
  agroup[0] = 1;
  agroup[1] = 1;
  elec.create(agroup);
  elec.R[0][0] = 0.00;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 0.0;
  elec.R[1][1] = 1.0;
  elec.R[1][2] = 0.0;

  EinsplineSetBuilder::PSetMap ptcl_map;
  ptcl_map.emplace(elec_ptr->getName(), std::move(elec_ptr));

  xmlNodePtr cur = NULL;
  EinsplineSetBuilder esb(elec, ptcl_map, c, cur);

  esb.SuperLattice       = 0.0;
  esb.SuperLattice(0, 0) = 1.0;
  esb.SuperLattice(1, 1) = 1.0;
  esb.SuperLattice(2, 2) = 1.0;

  REQUIRE(esb.CheckLattice());

  esb.SuperLattice(0, 0) = 1.1;
  REQUIRE_FALSE(esb.CheckLattice());
}

} // namespace qmcplusplus
