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
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"


#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Einspline SPO from HDF", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;


  elec_.setName("elec");
  elec_.create(2);
  elec_.R[0][0] = 0.00;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  // monoO
  /*
  elec_.Lattice.R(0,0) = 5.10509515;
  elec_.Lattice.R(0,1) = -3.23993545;
  elec_.Lattice.R(0,2) = 0.0;
  elec_.Lattice.R(1,0) = 5.10509515;
  elec_.Lattice.R(1,1) = 3.23993545;
  elec_.Lattice.R(1,2) = 0.0;
  elec_.Lattice.R(2,0) = -6.49690625;
  elec_.Lattice.R(2,1) = 0.0;
  elec_.Lattice.R(2,2) = 7.08268015;
 */

  // diamondC_1x1x1
  elec_.Lattice.R(0, 0) = 3.37316115;
  elec_.Lattice.R(0, 1) = 3.37316115;
  elec_.Lattice.R(0, 2) = 0.0;
  elec_.Lattice.R(1, 0) = 0.0;
  elec_.Lattice.R(1, 1) = 3.37316115;
  elec_.Lattice.R(1, 2) = 3.37316115;
  elec_.Lattice.R(2, 0) = 3.37316115;
  elec_.Lattice.R(2, 1) = 0.0;
  elec_.Lattice.R(2, 2) = 3.37316115;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

#ifdef ENABLE_SOA
  elec_.addTable(ions_, DT_SOA);
#else
  elec_.addTable(ions_, DT_AOS);
#endif
  elec_.resetGroups();
  elec_.update();


  TrialWaveFunction psi(c);
  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  SPOSet* spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != NULL);

#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // for vgl
  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
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
  SPOSet::ValueVector_t psiV(psiM[1], spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsiV(dpsiM[1], spo->getOrbitalSetSize());
  SPOSet::HessVector_t ddpsiV(spo->getOrbitalSetSize());
  spo->evaluate(elec_, 1, psiV, dpsiV, ddpsiV);

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

  SPOSet::HessMatrix_t hesspsiV(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GGGMatrix_t d3psiV(elec_.R.size(), spo->getOrbitalSetSize());
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
        SPOSet::ValueVector_t orbs(orbSize);
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

TEST_CASE("EinsplineSetBuilder CheckLattice", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet* elec = new ParticleSet;

  elec->setName("elec");
  std::vector<int> agroup(2);
  agroup[0] = 1;
  agroup[1] = 1;
  elec->create(agroup);
  elec->R[0][0] = 0.00;
  elec->R[0][1] = 0.0;
  elec->R[0][2] = 0.0;
  elec->R[1][0] = 0.0;
  elec->R[1][1] = 1.0;
  elec->R[1][2] = 0.0;

  elec->Lattice.R       = 0.0;
  elec->Lattice.R(0, 0) = 1.0;
  elec->Lattice.R(1, 1) = 1.0;
  elec->Lattice.R(2, 2) = 1.0;

  EinsplineSetBuilder::PtclPoolType ptcl_map;
  ptcl_map["e"] = elec;

  xmlNodePtr cur = NULL;
  EinsplineSetBuilder esb(*elec, ptcl_map, c, cur);

  esb.SuperLattice       = 0.0;
  esb.SuperLattice(0, 0) = 1.0;
  esb.SuperLattice(1, 1) = 1.0;
  esb.SuperLattice(2, 2) = 1.0;

  REQUIRE(esb.CheckLattice());

  esb.SuperLattice(0, 0) = 1.1;
  REQUIRE_FALSE(esb.CheckLattice());
}
} // namespace qmcplusplus
