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
#include "BsplineFactory/EinsplineSetBuilder.h"
#include "BsplineFactory/EinsplineSpinorSetBuilder.h"
#include <ResourceCollection.h>

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

  //diamondC_1x1x1
  // add save_coefs="yes" to create an HDF file of spline coefficients for the eval_bspline_spo.py script
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="8"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  // Test the case where the number of psi values is not the orbital set size
  // or the number of electrons/2
  const int psi_size = 3;
  SPOSet::ValueMatrix psiM(elec_.R.size(), psi_size);
  SPOSet::GradMatrix dpsiM(elec_.R.size(), psi_size);
  SPOSet::ValueMatrix d2psiM(elec_.R.size(), psi_size);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM, dpsiM, d2psiM);

  // value
  CHECK(std::real(psiM[0][0]) == Approx(-0.42546836868));
  CHECK(std::real(psiM[0][1]) == Approx(0.0));
  CHECK(std::real(psiM[1][0]) == Approx(-0.8886948824));
  CHECK(std::real(psiM[1][1]) == Approx(1.419412370359));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(-0.0000183403));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(0.1655139178));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0000193077));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(-1.3131694794));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(-1.1174004078));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.8462534547));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(1.3313053846));
  CHECK(std::real(d2psiM[1][1]) == Approx(-4.712583065));

  // for vgh
  SPOSet::ValueVector psiV(psiM[1], psi_size);
  SPOSet::GradVector dpsiV(dpsiM[1], psi_size);
  SPOSet::HessVector ddpsiV(psi_size);
  spo->evaluateVGH(elec_, 1, psiV, dpsiV, ddpsiV);

  // Catch default is 100*(float epsilson)
  double eps = 2000 * std::numeric_limits<float>::epsilon();
  //hess
  CHECK(std::real(ddpsiV[1](0, 0)) == Approx(-2.3160984034));
  CHECK(std::real(ddpsiV[1](0, 1)) == Approx(1.8089479397));
  CHECK(std::real(ddpsiV[1](0, 2)) == Approx(0.5608575749));
  CHECK(std::real(ddpsiV[1](1, 0)) == Approx(1.8089479397));
  CHECK(std::real(ddpsiV[1](1, 1)) == Approx(-0.07996207476).epsilon(eps));
  CHECK(std::real(ddpsiV[1](1, 2)) == Approx(0.5237969314));
  CHECK(std::real(ddpsiV[1](2, 0)) == Approx(0.5608575749));
  CHECK(std::real(ddpsiV[1](2, 1)) == Approx(0.5237969314));
  CHECK(std::real(ddpsiV[1](2, 2)) == Approx(-2.316497764));

  SPOSet::HessMatrix hesspsiV(elec_.R.size(), psi_size);
  SPOSet::GGGMatrix d3psiV(elec_.R.size(), psi_size);
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
  CHECK(std::real(d3psiV(1,0)[0][0] ) == Approx(0.0463371276));
  CHECK(std::real(d3psiV(1,0)[0][1] ) == Approx(1.1755813360));
  CHECK(std::real(d3psiV(1,0)[0][2] ) == Approx(0.0660155713));
  CHECK(std::real(d3psiV(1,0)[0][4] ) == Approx(0.0414704382));
  CHECK(std::real(d3psiV(1,0)[0][5] ) == Approx(-0.5167412758));
  CHECK(std::real(d3psiV(1,0)[0][8] ) == Approx(0.0659536421));
  CHECK(std::real(d3psiV(1,0)[1][4] ) == Approx(-4.8771157264));
  CHECK(std::real(d3psiV(1,0)[1][5] ) == Approx(0.0415326356));
  CHECK(std::real(d3psiV(1,0)[1][8] ) == Approx(1.1755810976));
  CHECK(std::real(d3psiV(1,0)[2][8] ) == Approx(0.0463993549));

  CHECK(std::real(d3psiV(1,1)[0][0] ) == Approx(6.7155771255));
  CHECK(std::real(d3psiV(1,1)[0][1] ) == Approx(5.5450510978));
  CHECK(std::real(d3psiV(1,1)[0][2] ) == Approx(0.9829711914));
  CHECK(std::real(d3psiV(1,1)[0][4] ) == Approx(-3.1704092025));
  CHECK(std::real(d3psiV(1,1)[0][5] ) == Approx(-1.9537661075));
  CHECK(std::real(d3psiV(1,1)[0][8] ) == Approx(1.9305641651));
  CHECK(std::real(d3psiV(1,1)[1][4] ) == Approx(3.6051378250));
  CHECK(std::real(d3psiV(1,1)[1][5] ) == Approx(-0.7382576465));
  CHECK(std::real(d3psiV(1,1)[1][8] ) == Approx(5.5741839408));
  CHECK(std::real(d3psiV(1,1)[2][8] ) == Approx(3.1312348842));
#endif

  // Test the batched interface
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

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

  const int ne = elec_.R.size();
  SPOSet::ValueMatrix psi(ne, psi_size);
  SPOSet::GradMatrix dpsi(ne, psi_size);
  SPOSet::ValueMatrix d2psi(ne, psi_size);
  SPOSet::ValueMatrix psi_2(ne, psi_size);
  SPOSet::GradMatrix dpsi_2(ne, psi_size);
  SPOSet::ValueMatrix d2psi_2(ne, psi_size);

  RefVector<SPOSet::ValueMatrix> psi_v_list;
  RefVector<SPOSet::GradMatrix> dpsi_v_list;
  RefVector<SPOSet::ValueMatrix> d2psi_v_list;

  psi_v_list.push_back(psi);
  psi_v_list.push_back(psi_2);
  dpsi_v_list.push_back(dpsi);
  dpsi_v_list.push_back(dpsi_2);
  d2psi_v_list.push_back(d2psi);
  d2psi_v_list.push_back(d2psi_2);
  spo->mw_evaluate_notranspose(spo_list, p_list, 0, elec_.R.size(), psi_v_list, dpsi_v_list, d2psi_v_list);

  CHECK(std::real(psi_v_list[0].get()[0][0]) == Approx(-0.42546836868));
  CHECK(std::real(psi_v_list[0].get()[1][0]) == Approx(-0.8886948824));
  // Second particle set had particle positions flipped
  CHECK(std::real(psi_v_list[1].get()[0][0]) == Approx(-0.8886948824));
  CHECK(std::real(psi_v_list[1].get()[1][0]) == Approx(-0.42546836868));

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
        elec_.R[0] = {x, y, z};
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

TEST_CASE("Einspline SPO from HDF diamond_2x1x1 5 electrons", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  // diamondC_2x1x1
  lattice.R = {6.7463223, 6.7463223, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({4});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  ions_.R[2] = {3.37316115, 3.37316115, 0.0};
  ions_.R[3] = {5.05974173, 5.05974173, 1.68658058};


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({5});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};
  elec_.R[2] = {0.0, 1.1, 0.0};
  elec_.R[3] = {0.0, 1.2, 0.0};
  elec_.R[4] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_2x1x1
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_2x1x1.pwscf.h5" tilematrix="2 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="5"/>
</tmp>
)";

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

  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  CHECK(std::real(psiM[1][0]) == Approx(0.9008999467));
  CHECK(std::real(psiM[1][1]) == Approx(1.2383049726));
  // grad
  CHECK(std::real(dpsiM[1][0][0]) == Approx(0.0025820041));
  CHECK(std::real(dpsiM[1][0][1]) == Approx(-0.1880052537));
  CHECK(std::real(dpsiM[1][0][2]) == Approx(-0.0025404284));
  CHECK(std::real(dpsiM[1][1][0]) == Approx(0.1069662273));
  CHECK(std::real(dpsiM[1][1][1]) == Approx(-0.4364597797));
  CHECK(std::real(dpsiM[1][1][2]) == Approx(-0.106951952));
  // lapl
  CHECK(std::real(d2psiM[1][0]) == Approx(-1.3757134676));
  CHECK(std::real(d2psiM[1][1]) == Approx(-2.4803137779));

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  CHECK(std::imag(psiM[1][0]) == Approx(0.9008999467));
  CHECK(std::imag(psiM[1][1]) == Approx(1.2383049726));
  // grad
  CHECK(std::imag(dpsiM[1][0][0]) == Approx(0.0025820041));
  CHECK(std::imag(dpsiM[1][0][1]) == Approx(-0.1880052537));
  CHECK(std::imag(dpsiM[1][0][2]) == Approx(-0.0025404284));
  CHECK(std::imag(dpsiM[1][1][0]) == Approx(0.1069453433));
  CHECK(std::imag(dpsiM[1][1][1]) == Approx(-0.43649593));
  CHECK(std::imag(dpsiM[1][1][2]) == Approx(-0.1069145575));
  // lapl
  CHECK(std::imag(d2psiM[1][0]) == Approx(-1.3757134676));
  CHECK(std::imag(d2psiM[1][1]) == Approx(-2.4919104576));
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

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec_.createResource(pset_res);
  spo->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

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
  // real part
  // due to the different ordering of bands skip the tests on CUDA+Real builds
  // checking evaluations, reference values are not independently generated.
  // value
  CHECK(std::real(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  CHECK(std::real(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  CHECK(std::real(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  CHECK(std::real(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  CHECK(std::real(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  CHECK(std::real(dpsi_v_list[1].get()[1][0]) == Approx(0.1069662273));
  CHECK(std::real(dpsi_v_list[1].get()[1][1]) == Approx(-0.4364597797));
  CHECK(std::real(dpsi_v_list[1].get()[1][2]) == Approx(-0.106951952));
  // lapl
  CHECK(std::real(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  CHECK(std::real(d2psi_v_list[1].get()[1]) == Approx(-2.4803137779));

#if defined(QMC_COMPLEX)
  // imaginary part
  // value
  CHECK(std::imag(psi_v_list[1].get()[0]) == Approx(0.9008999467));
  CHECK(std::imag(psi_v_list[1].get()[1]) == Approx(1.2383049726));
  // grad
  CHECK(std::imag(dpsi_v_list[1].get()[0][0]) == Approx(0.0025820041));
  CHECK(std::imag(dpsi_v_list[1].get()[0][1]) == Approx(-0.1880052537));
  CHECK(std::imag(dpsi_v_list[1].get()[0][2]) == Approx(-0.0025404284));
  CHECK(std::imag(dpsi_v_list[1].get()[1][0]) == Approx(0.1069453433));
  CHECK(std::imag(dpsi_v_list[1].get()[1][1]) == Approx(-0.43649593));
  CHECK(std::imag(dpsi_v_list[1].get()[1][2]) == Approx(-0.1069145575));
  // lapl
  CHECK(std::imag(d2psi_v_list[1].get()[0]) == Approx(-1.3757134676));
  CHECK(std::imag(d2psi_v_list[1].get()[1]) == Approx(-2.4919104576));
#endif

  const size_t nw = 2;
  std::vector<SPOSet::ValueType> ratio_v(nw);
  std::vector<SPOSet::GradType> grads_v(nw);

  Vector<SPOSet::ValueType, OffloadPinnedAllocator<SPOSet::ValueType>> inv_row(5);
  inv_row = {0.1, 0.2, 0.3, 0.4, 0.5};
  inv_row.updateTo();

  std::vector<const SPOSet::ValueType*> inv_row_ptr(nw, inv_row.device_data());

  SPOSet::OffloadMWVGLArray phi_vgl_v;
  phi_vgl_v.resize(QMCTraits::DIM_VGL, nw, 5);
  spo->mw_evaluateVGLandDetRatioGrads(spo_list, p_list, 0, inv_row_ptr, phi_vgl_v, ratio_v, grads_v);
#if !defined(QMC_COMPLEX)
  CHECK(std::real(ratio_v[0]) == Approx(0.2365307168));
  CHECK(std::real(grads_v[0][0]) == Approx(-5.4095164399));
  CHECK(std::real(grads_v[0][1]) == Approx(14.37990087));
  CHECK(std::real(grads_v[0][2]) == Approx(16.9374788259));
  CHECK(std::real(ratio_v[1]) == Approx(1.0560744941));
  CHECK(std::real(grads_v[1][0]) == Approx(-0.0863436466));
  CHECK(std::real(grads_v[1][1]) == Approx(-0.7499371447));
  CHECK(std::real(grads_v[1][2]) == Approx(0.8570534314));
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
  elec.R[0] = {0.0, 0.0, 0.0};
  elec.R[1] = {0.0, 1.0, 0.0};

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
