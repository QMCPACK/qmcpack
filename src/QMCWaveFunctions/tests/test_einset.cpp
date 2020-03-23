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
  elec_.R[0][0] = 0.0;
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
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;

  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  //diamondC_1x1x1
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"4\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  std::unique_ptr<SPOSet> spo(einSet.createSPOSetFromXML(ein1));
  REQUIRE(spo);

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

TEST_CASE("Einspline SPO from HDF diamond_2x1x1", "[wavefunction]")
{
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(4);
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
  elec_.create(2);
  elec_.R[0][0] = 0.0;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  // diamondC_2x1x1
  elec_.Lattice.R(0, 0) = 6.7463223;
  elec_.Lattice.R(0, 1) = 6.7463223;
  elec_.Lattice.R(0, 2) = 0.0;
  elec_.Lattice.R(1, 0) = 0.0;
  elec_.Lattice.R(1, 1) = 3.37316115;
  elec_.Lattice.R(1, 2) = 3.37316115;
  elec_.Lattice.R(2, 0) = 3.37316115;
  elec_.Lattice.R(2, 1) = 0.0;
  elec_.Lattice.R(2, 2) = 3.37316115;

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;

  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

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
  std::unique_ptr<SPOSet> spo(einSet.createSPOSetFromXML(ein1));
  REQUIRE(spo);

  // for vgl
  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
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
  std::vector<ParticleSet*> P_list;
  P_list.push_back(&elec_);
  P_list.push_back(&elec_2);

  std::unique_ptr<SPOSet> spo_2(spo->makeClone());
  std::vector<SPOSet*> spo_list;
  spo_list.push_back(spo.get());
  spo_list.push_back(spo_2.get());

  SPOSet::ValueVector_t psi(spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsi(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t d2psi(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t psi_2(spo->getOrbitalSetSize());
  SPOSet::GradVector_t dpsi_2(spo->getOrbitalSetSize());
  SPOSet::ValueVector_t d2psi_2(spo->getOrbitalSetSize());

  RefVector<SPOSet::ValueVector_t> psi_v_list;
  RefVector<SPOSet::GradVector_t> dpsi_v_list;
  RefVector<SPOSet::ValueVector_t> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);

  spo->mw_evaluateValue(spo_list, P_list, 1, psi_v_list);
#if !defined(QMC_CUDA) || defined(QMC_COMPLEX)
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  REQUIRE(std::real(psi_v_list[0].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::real(psi_v_list[0].get()[1]) == Approx(1.2383049726));
#endif

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  REQUIRE(std::imag(psi_v_list[0].get()[0]) == Approx(0.9008999467));
  REQUIRE(std::imag(psi_v_list[0].get()[1]) == Approx(1.2383049726));
#endif

  spo->mw_evaluateVGL(spo_list, P_list, 0, psi_v_list, dpsi_v_list, d2psi_v_list);
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

//Now we test the spinor set with Einspline orbitals from HDF.
#ifdef QMC_COMPLEX
TEST_CASE("Einspline SpinorSet from HDF", "[wavefunction]")
{
  app_log()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log()<<"!!!!!  Einspline SpinorSet from HDF   !!!!!\n";
  app_log()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
 
  using ValueType=SPOSet::ValueType;
  using RealType=SPOSet::RealType;
  Communicate* c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);

  ions_.R[0][0] = 0.00000000 ; 
  ions_.R[0][1] = 0.00000000 ;    
  ions_.R[0][2] = 1.08659253 ;  
  ions_.R[1][0] = 0.00000000 ;   
  ions_.R[1][1] = 0.00000000 ;  
  ions_.R[1][2] =-1.08659253 ;

  elec_.setName("elec");
  elec_.create(3);
  elec_.R[0][0] = 0.1;
  elec_.R[0][1] = -0.3;
  elec_.R[0][2] = 1.0;
  elec_.R[1][0] = -0.1;
  elec_.R[1][1] = 0.3;
  elec_.R[1][2] = 1.0;
  elec_.R[2][0] = 0.1;
  elec_.R[2][1] = 0.2;
  elec_.R[2][2] = 0.3;

  elec_.spins[0] = 0.0 ;
  elec_.spins[1] = 0.2 ;
  elec_.spins[2] = 0.4 ;

  // O2 test example from pwscf non-collinear calculation.
  elec_.Lattice.R(0, 0) =  5.10509515 ;
  elec_.Lattice.R(0, 1) = -3.23993545 ;
  elec_.Lattice.R(0, 2) =  0.00000000 ;
  elec_.Lattice.R(1, 0) = 5.10509515 ;
  elec_.Lattice.R(1, 1) = 3.23993545 ;
  elec_.Lattice.R(1, 2) = 0.00000000 ;
  elec_.Lattice.R(2, 0) = -6.49690625 ;
  elec_.Lattice.R(2, 1) =  0.00000000 ;
  elec_.Lattice.R(2, 2) =  7.08268015 ; 

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);


  const char* particles = "<tmp> \
   <sposet_builder name=\"A\" type=\"spinorbspline\" href=\"o2_45deg_spins.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" size=\"3\" precision=\"float\"> \
     <sposet name=\"myspo\" size=\"3\"> \
       <occupation mode=\"ground\"/> \
     </sposet> \
   </sposet_builder> \
   </tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSpinorSetBuilder einSet(elec_,ptcl.getPool(), c, ein1);
  std::unique_ptr<SPOSet> spo(einSet.createSPOSetFromXML(ein1));
  REQUIRE(spo);

  SPOSet::ValueMatrix_t psiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM(elec_.R.size(), spo->getOrbitalSetSize()); //spin gradient
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
 
  //These are the reference values computed from a spin-polarized calculation, 
  //with the assumption that the coefficients for phi^\uparrow 
  SPOSet::ValueMatrix_t psiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_ref(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t dspsiM_ref(elec_.R.size(), spo->getOrbitalSetSize());

  SPOSet::GradMatrix_t dpsiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::GradMatrix_t dpsiM_ref(elec_.R.size(), spo->getOrbitalSetSize());

  SPOSet::ValueMatrix_t d2psiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t d2psiM_ref(elec_.R.size(), spo->getOrbitalSetSize());


  //These reference values were generated as follows:
  // 1.) Non-Collinear O2 calculation in PBC's performed using Quantum Espresso.  
  // 2.) Spinor wavefunction converted to HDF5 using convertpw4qmc tool.  Mainly, this places the up channel in the spin_0 slot, and spin down in spin_1.  
  // 3.) The HDF5 metadata was hacked by hand to correspond to a fictional but consistent spin-polarized set of orbitals.
  // 4.) A spin polarized QMCPACK run was done using the electron and ion configuration specified in this block.  Orbital values, gradients, and laplacians were calculated
  //     for both "spin up" and "spin down" orbitals. This is where the psiM_up(down), d2psiM_up(down) values come from.   
  // 5.) By hand, the reference values, gradients, and laplacians are calculated by using the formula for a spinor e^is phi_up + e^{-is} phi_down.
  // 6.) These are compared against the integrated initialization/parsing/evaluation of the Einspline Spinor object.  
  
  //Reference values for spin up component.
  psiM_up[0][0]=ValueType(2.8696985245e+00,-2.8696982861e+00); psiM_up[0][1]=ValueType(1.1698637009e+00,-1.1698638201e+00); psiM_up[0][2]=ValueType(-2.6149117947e+00,2.6149117947e+00);
  psiM_up[1][0]=ValueType(2.8670933247e+00,-2.8670933247e+00); psiM_up[1][1]=ValueType(1.1687355042e+00,-1.1687356234e+00); psiM_up[1][2]=ValueType(-2.6131081581e+00,2.6131081581e+00);
  psiM_up[2][0]=ValueType(4.4833350182e+00,-4.4833350182e+00); psiM_up[2][1]=ValueType(1.8927993774e+00,-1.8927993774e+00); psiM_up[2][2]=ValueType(-8.3977413177e-01,8.3977431059e-01);

  //Reference values for spin down component.
  psiM_down[0][0]=ValueType(1.1886650324e+00,-1.1886655092e+00); psiM_down[0][1]=ValueType(-2.8243079185e+00,2.8243076801e+00); psiM_down[0][2]=ValueType(-1.0831292868e+00,1.0831292868e+00);
  psiM_down[1][0]=ValueType(1.1875861883e+00,-1.1875866652e+00); psiM_down[1][1]=ValueType(-2.8215842247e+00,2.8215837479e+00); psiM_down[1][2]=ValueType(-1.0823822021e+00,1.0823823214e+00);
  psiM_down[2][0]=ValueType(1.8570541143e+00,-1.8570543528e+00); psiM_down[2][1]=ValueType(-4.5696320534e+00,4.5696320534e+00); psiM_down[2][2]=ValueType(-3.4784498811e-01,3.4784474969e-01);

  //And the laplacians...
  d2psiM_up[0][0]=ValueType(-6.1587309837e+00,6.1587429047e+00); d2psiM_up[0][1]=ValueType(-2.4736759663e+00,2.4736781120e+00); d2psiM_up[0][2]=ValueType(2.1381640434e-01,-2.1381306648e-01);
  d2psiM_up[1][0]=ValueType(-5.0561609268e+00,5.0561575890e+00); d2psiM_up[1][1]=ValueType(-2.0328726768e+00,2.0328762531e+00); d2psiM_up[1][2]=ValueType(-7.4090242386e-01,7.4090546370e-01);
  d2psiM_up[2][0]=ValueType(-1.8970542908e+01,1.8970539093e+01); d2psiM_up[2][1]=ValueType(-8.2134075165e+00,8.2134037018e+00); d2psiM_up[2][2]=ValueType(1.0161912441e+00,-1.0161914825e+00);

  d2psiM_down[0][0]=ValueType(-2.5510206223e+00,2.5510258675e+00); d2psiM_down[0][1]=ValueType(5.9720201492e+00,-5.9720129967e+00); d2psiM_down[0][2]=ValueType(8.8568925858e-02,-8.8571548462e-02);
  d2psiM_down[1][0]=ValueType(-2.0943276882e+00,2.0943336487e+00); d2psiM_down[1][1]=ValueType(4.9078116417e+00,-4.9078197479e+00); d2psiM_down[1][2]=ValueType(-3.0689623952e-01,3.0689093471e-01);
  d2psiM_down[2][0]=ValueType(-7.8578405380e+00,7.8578381538e+00); d2psiM_down[2][1]=ValueType(1.9828968048e+01,-1.9828992844e+01); d2psiM_down[2][2]=ValueType(4.2092007399e-01,-4.2091816664e-01);

  //And now a looooot of gradient info.
 ////////SPIN UP//////////
  dpsiM_up[0][0][0]=ValueType(-1.7161563039e-01,1.7161482573e-01);   dpsiM_up[0][0][1]=ValueType(5.6693041325e-01,-5.6692999601e-01);   dpsiM_up[0][0][2]=ValueType(-4.5538558960e+00,4.5538554192e+00); 
  dpsiM_up[0][1][0]=ValueType(-7.4953302741e-02,7.4952393770e-02);   dpsiM_up[0][1][1]=ValueType(2.4608184397e-01,-2.4608163536e-01);   dpsiM_up[0][1][2]=ValueType(-1.9720511436e+00,1.9720509052e+00); 
  dpsiM_up[0][2][0]=ValueType(-4.2384520173e-02,4.2384237051e-02);   dpsiM_up[0][2][1]=ValueType(1.1735939980e-01,-1.1735984683e-01);   dpsiM_up[0][2][2]=ValueType(-3.1189033985e+00,3.1189031601e+00); 
 
  dpsiM_up[1][0][0]=ValueType(1.9333077967e-01,-1.9333113730e-01);   dpsiM_up[1][0][1]=ValueType(-5.7470333576e-01,5.7470202446e-01);   dpsiM_up[1][0][2]=ValueType(-4.5568108559e+00,4.5568113327e+00); 
  dpsiM_up[1][1][0]=ValueType(8.4540992975e-02,-8.4540143609e-02);   dpsiM_up[1][1][1]=ValueType(-2.4946013093e-01,2.4946044385e-01);   dpsiM_up[1][1][2]=ValueType(-1.9727530479e+00,1.9727528095e+00); 
  dpsiM_up[1][2][0]=ValueType(3.1103719026e-02,-3.1103719026e-02);   dpsiM_up[1][2][1]=ValueType(-1.2540178001e-01,1.2540178001e-01);   dpsiM_up[1][2][2]=ValueType(-3.1043677330e+00,3.1043677330e+00); 
 
  dpsiM_up[2][0][0]=ValueType(-8.8733488321e-01,8.8733488321e-01);   dpsiM_up[2][0][1]=ValueType(-1.7726477385e+00,1.7726477385e+00);   dpsiM_up[2][0][2]=ValueType(7.3728728294e-01,-7.3728692532e-01); 
  dpsiM_up[2][1][0]=ValueType(-3.8018247485e-01,3.8018330932e-01);   dpsiM_up[2][1][1]=ValueType(-7.5880718231e-01,7.5880759954e-01);   dpsiM_up[2][1][2]=ValueType(2.7537062764e-01,-2.7537041903e-01); 
  dpsiM_up[2][2][0]=ValueType(-9.5389984548e-02,9.5390148461e-02);   dpsiM_up[2][2][1]=ValueType(-1.8467208743e-01,1.8467210233e-01);   dpsiM_up[2][2][2]=ValueType(-2.4704084396e+00,2.4704084396e+00); 
 
 ////////SPIN DOWN//////////
  dpsiM_down[0][0][0]=ValueType(-7.1084961295e-02,7.1085616946e-02);   dpsiM_down[0][0][1]=ValueType(2.3483029008e-01,-2.3482969403e-01);   dpsiM_down[0][0][2]=ValueType(-1.8862648010e+00,1.8862643242e+00); 
  dpsiM_down[0][1][0]=ValueType(1.8095153570e-01,-1.8095159531e-01);   dpsiM_down[0][1][1]=ValueType(-5.9409534931e-01,5.9409546852e-01);   dpsiM_down[0][1][2]=ValueType(4.7609643936e+00,-4.7609624863e+00); 
  dpsiM_down[0][2][0]=ValueType(-1.7556600273e-02,1.7556769773e-02);   dpsiM_down[0][2][1]=ValueType(4.8611730337e-02,-4.8612065613e-02);   dpsiM_down[0][2][2]=ValueType(-1.2918885946e+00,1.2918891907e+00); 
 
  dpsiM_down[1][0][0]=ValueType(8.0079451203e-02,-8.0079004169e-02);   dpsiM_down[1][0][1]=ValueType(-2.3804906011e-01,2.3804855347e-01);   dpsiM_down[1][0][2]=ValueType(-1.8874882460e+00,1.8874886036e+00); 
  dpsiM_down[1][1][0]=ValueType(-2.0409825444e-01,2.0409949124e-01);   dpsiM_down[1][1][1]=ValueType(6.0225284100e-01,-6.0225236416e-01);   dpsiM_down[1][1][2]=ValueType(4.7626581192e+00,-4.7626576424e+00); 
  dpsiM_down[1][2][0]=ValueType(1.2884057127e-02,-1.2884397991e-02);   dpsiM_down[1][2][1]=ValueType(-5.1943197846e-02,5.1943652332e-02);   dpsiM_down[1][2][2]=ValueType(-1.2858681679e+00,1.2858685255e+00); 
 
  dpsiM_down[2][0][0]=ValueType(-3.6754500866e-01,3.6754477024e-01);   dpsiM_down[2][0][1]=ValueType(-7.3425340652e-01,7.3425388336e-01);   dpsiM_down[2][0][2]=ValueType(3.0539327860e-01,-3.0539402366e-01); 
  dpsiM_down[2][1][0]=ValueType(9.1784656048e-01,-9.1784602404e-01);   dpsiM_down[2][1][1]=ValueType(1.8319253922e+00,-1.8319247961e+00);   dpsiM_down[2][1][2]=ValueType(-6.6480386257e-01,6.6480308771e-01); 
  dpsiM_down[2][2][0]=ValueType(-3.9511863142e-02,3.9511814713e-02);   dpsiM_down[2][2][1]=ValueType(-7.6493337750e-02,7.6493576169e-02);   dpsiM_down[2][2][2]=ValueType(-1.0232743025e+00,1.0232743025e+00); 


  for(unsigned int iat=0; iat<3; iat++)
  {
    RealType s = elec_.spins[iat];
    RealType coss(0.0),sins(0.0);

    coss=std::cos(s);
    sins=std::sin(s);

    ValueType eis(coss,sins);
    ValueType emis(coss,-sins);
    ValueType eye(0,1.0);
    //Using the reference values for the up and down channels invdividually, we build the total reference spinor value
    //consistent with the current spin value of particle iat.
    psiM_ref[iat][0]=eis*psiM_up[iat][0]+emis*psiM_down[iat][0];
    psiM_ref[iat][1]=eis*psiM_up[iat][1]+emis*psiM_down[iat][1];
    psiM_ref[iat][2]=eis*psiM_up[iat][2]+emis*psiM_down[iat][2];
  
    dspsiM_ref[iat][0] = eye*(eis*psiM_up[iat][0] - emis*psiM_down[iat][0]);
    dspsiM_ref[iat][1] = eye*(eis*psiM_up[iat][1] - emis*psiM_down[iat][1]);
    dspsiM_ref[iat][2] = eye*(eis*psiM_up[iat][2] - emis*psiM_down[iat][2]);
  
    dpsiM_ref[iat][0]=eis*dpsiM_up[iat][0]+emis*dpsiM_down[iat][0];
    dpsiM_ref[iat][1]=eis*dpsiM_up[iat][1]+emis*dpsiM_down[iat][1];
    dpsiM_ref[iat][2]=eis*dpsiM_up[iat][2]+emis*dpsiM_down[iat][2];
    
    d2psiM_ref[iat][0]=eis*d2psiM_up[iat][0]+emis*d2psiM_down[iat][0];
    d2psiM_ref[iat][1]=eis*d2psiM_up[iat][1]+emis*d2psiM_down[iat][1];
    d2psiM_ref[iat][2]=eis*d2psiM_up[iat][2]+emis*d2psiM_down[iat][2];

  }

  //OK.  Going to test evaluate_notranspose with spin.
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);
 
  for(unsigned int iat=0; iat<3; iat++)
  {
    
    REQUIRE( psiM[iat][0] == ComplexApprox( psiM_ref[iat][0] ));
    REQUIRE( psiM[iat][1] == ComplexApprox( psiM_ref[iat][1] ));
    REQUIRE( psiM[iat][2] == ComplexApprox( psiM_ref[iat][2] ));

    REQUIRE( dpsiM[iat][0][0] == ComplexApprox( dpsiM_ref[iat][0][0] ));
    REQUIRE( dpsiM[iat][0][1] == ComplexApprox( dpsiM_ref[iat][0][1] ));
    REQUIRE( dpsiM[iat][0][2] == ComplexApprox( dpsiM_ref[iat][0][2] ));

    REQUIRE( dpsiM[iat][1][0] == ComplexApprox( dpsiM_ref[iat][1][0] ));
    REQUIRE( dpsiM[iat][1][1] == ComplexApprox( dpsiM_ref[iat][1][1] ));
    REQUIRE( dpsiM[iat][1][2] == ComplexApprox( dpsiM_ref[iat][1][2] ));
    
    REQUIRE( dpsiM[iat][2][0] == ComplexApprox( dpsiM_ref[iat][2][0] ));
    REQUIRE( dpsiM[iat][2][1] == ComplexApprox( dpsiM_ref[iat][2][1] ));
    REQUIRE( dpsiM[iat][2][2] == ComplexApprox( dpsiM_ref[iat][2][2] ));



    REQUIRE( d2psiM[iat][0] == ComplexApprox( d2psiM_ref[iat][0] ));
    REQUIRE( d2psiM[iat][1] == ComplexApprox( d2psiM_ref[iat][1] ));
    REQUIRE( d2psiM[iat][2] == ComplexApprox( d2psiM_ref[iat][2] ));
  }

  //Now we're going to test evaluateValue and evaluateVGL:

  int OrbitalSetSize = spo->getOrbitalSetSize();  
  //temporary arrays for holding the values of the up and down channels respectively.
  SPOSet::ValueVector_t psi_work;

  //temporary arrays for holding the gradients of the up and down channels respectively.
  SPOSet::GradVector_t dpsi_work;

  //temporary arrays for holding the laplacians of the up and down channels respectively.
  SPOSet::ValueVector_t d2psi_work;

   
  psi_work.resize(OrbitalSetSize);

  dpsi_work.resize(OrbitalSetSize);

  d2psi_work.resize(OrbitalSetSize);

  //We worked hard to generate nice reference data above.  Let's generate a test for evaluateV 
  //and evaluateVGL by perturbing the electronic configuration by dR, and then make
  //single particle moves that bring it back to our original R reference values.

  //Our perturbation vector. 
  ParticleSet::ParticlePos_t dR;
  dR.resize(3);
  //Values chosen based on divine inspiration.  Not important.
  dR[0][0] = 0.1; dR[0][1] = 0.2; dR[0][2] = 0.1;
  dR[1][0] = -0.1; dR[1][1] = -0.05 ; dR[1][2] = 0.05;
  dR[2][0] = -0.1; dR[2][1] =0.0 ; dR[2][2] =0.0 ;

  //The new R of our perturbed particle set.
  ParticleSet::ParticlePos_t Rnew;
  Rnew.resize(3);
  Rnew = elec_.R+dR;
  elec_.R=Rnew;
  elec_.update(); 
 
  //Now we test evaluateValue()
  for(unsigned int iat=0; iat<3; iat++)
  {
    psi_work=0.0;
    elec_.makeMove(iat,-dR[iat],false);
    spo->evaluateValue(elec_,iat,psi_work);
  
    REQUIRE( psi_work[0] == ComplexApprox( psiM_ref[iat][0] ));
    REQUIRE( psi_work[1] == ComplexApprox( psiM_ref[iat][1] ));
    REQUIRE( psi_work[2] == ComplexApprox( psiM_ref[iat][2] ));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateVGL()
  for(unsigned int iat=0; iat<3; iat++)
  {
    psi_work=0.0;
    dpsi_work=0.0;
    d2psi_work=0.0;

    elec_.makeMove(iat,-dR[iat],false);
    spo->evaluateVGL(elec_,iat,psi_work,dpsi_work,d2psi_work);
  
    REQUIRE( psi_work[0] == ComplexApprox( psiM_ref[iat][0] ));
    REQUIRE( psi_work[1] == ComplexApprox( psiM_ref[iat][1] ));
    REQUIRE( psi_work[2] == ComplexApprox( psiM_ref[iat][2] ));

    REQUIRE( dpsi_work[0][0] == ComplexApprox( dpsiM_ref[iat][0][0] ));
    REQUIRE( dpsi_work[0][1] == ComplexApprox( dpsiM_ref[iat][0][1] ));
    REQUIRE( dpsi_work[0][2] == ComplexApprox( dpsiM_ref[iat][0][2] ));

    REQUIRE( dpsi_work[1][0] == ComplexApprox( dpsiM_ref[iat][1][0] ));
    REQUIRE( dpsi_work[1][1] == ComplexApprox( dpsiM_ref[iat][1][1] ));
    REQUIRE( dpsi_work[1][2] == ComplexApprox( dpsiM_ref[iat][1][2] ));

    REQUIRE( dpsi_work[2][0] == ComplexApprox( dpsiM_ref[iat][2][0] ));
    REQUIRE( dpsi_work[2][1] == ComplexApprox( dpsiM_ref[iat][2][1] ));
    REQUIRE( dpsi_work[2][2] == ComplexApprox( dpsiM_ref[iat][2][2] ));

    REQUIRE( d2psi_work[0] == ComplexApprox( d2psiM_ref[iat][0] ));
    REQUIRE( d2psi_work[1] == ComplexApprox( d2psiM_ref[iat][1] ));
    REQUIRE( d2psi_work[2] == ComplexApprox( d2psiM_ref[iat][2] ));
    elec_.rejectMove(iat);
  }

  //Now we test evaluateSpin:
  SPOSet::ValueVector_t dspsi_work;
  dspsi_work.resize(OrbitalSetSize);

  for(unsigned int iat=0; iat<3; iat++)
  {
    psi_work=0.0;
    dspsi_work=0.0;
    
    elec_.makeMove(iat,-dR[iat],false);
    spo->evaluate_spin(elec_,iat,psi_work,dspsi_work);
    
    REQUIRE( psi_work[0] == ComplexApprox( psiM_ref[iat][0] ));
    REQUIRE( psi_work[1] == ComplexApprox( psiM_ref[iat][1] ));
    REQUIRE( psi_work[2] == ComplexApprox( psiM_ref[iat][2] ));

    REQUIRE( dspsi_work[0] == ComplexApprox( dspsiM_ref[iat][0] ));
    REQUIRE( dspsi_work[1] == ComplexApprox( dspsiM_ref[iat][1] ));
    REQUIRE( dspsi_work[2] == ComplexApprox( dspsiM_ref[iat][2] ));

    elec_.rejectMove(iat);
  }

}
#endif //QMC_COMPLEX


} // namespace qmcplusplus
