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
#include "QMCApp/ParticleSetPool.h"
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
  SPOSet::ValueMatrix_t d2psiM(elec_.R.size(), spo->getOrbitalSetSize());
 
  //These are the reference values computed from a spin-polarized calculation, 
  //with the assumption that the coefficients for phi^\uparrow 
  SPOSet::ValueMatrix_t psiM_up(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_down(elec_.R.size(), spo->getOrbitalSetSize());
  SPOSet::ValueMatrix_t psiM_ref(elec_.R.size(), spo->getOrbitalSetSize());

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

  for(unsigned int iat=0; iat<3; iat++)
  {
    RealType s = elec_.spins[iat];
    RealType coss(0.0),sins(0.0);

    coss=std::cos(s);
    sins=std::sin(s);

    ValueType eis(coss,sins);
    ValueType emis(coss,-sins);
    //Using the reference values for the up and down channels invdividually, we build the total reference spinor value
    //consistent with the current spin value of particle iat.
    psiM_ref[iat][0]=eis*psiM_up[iat][0]+emis*psiM_down[iat][0];
    psiM_ref[iat][1]=eis*psiM_up[iat][1]+emis*psiM_down[iat][1];
    psiM_ref[iat][2]=eis*psiM_up[iat][2]+emis*psiM_down[iat][2];

    d2psiM_ref[iat][0]=eis*d2psiM_up[iat][0]+emis*d2psiM_down[iat][0];
    d2psiM_ref[iat][1]=eis*d2psiM_up[iat][1]+emis*d2psiM_down[iat][1];
    d2psiM_ref[iat][2]=eis*d2psiM_up[iat][2]+emis*d2psiM_down[iat][2];

  }

  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);
 
  for(unsigned int iat=0; iat<3; iat++)
  {
    
    REQUIRE( psiM[iat][0] == ComplexApprox( psiM_ref[iat][0] ));
    REQUIRE( psiM[iat][1] == ComplexApprox( psiM_ref[iat][1] ));
    REQUIRE( psiM[iat][2] == ComplexApprox( psiM_ref[iat][2] ));

    REQUIRE( d2psiM[iat][0] == ComplexApprox( d2psiM_ref[iat][0] ));
    REQUIRE( d2psiM[iat][1] == ComplexApprox( d2psiM_ref[iat][1] ));
    REQUIRE( d2psiM[iat][2] == ComplexApprox( d2psiM_ref[iat][2] ));
  }

}
#endif //QMC_COMPLEX


} // namespace qmcplusplus
