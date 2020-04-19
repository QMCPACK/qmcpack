//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"

namespace qmcplusplus
{

using DiracDet = DiracDeterminantBatched<MatrixUpdateOMP<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;

using LogValueType = TrialWaveFunction::LogValueType;
using PsiValueType = TrialWaveFunction::PsiValueType;
using GradType = TrialWaveFunction::GradType;

TEST_CASE("TrialWaveFunction_diamondC_2x1x1", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

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
  ions_.update();


  elec_.setName("elec");
  std::vector<int> ud(2);
  ud[0] = ud[1] = 2;
  elec_.create(ud);
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 1.0};
  elec_.R[2] = {1.0, 1.0, 0.0};
  elec_.R[3] = {1.0, 0.0, 1.0};

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
  elec_.Lattice.BoxBConds = {1, 1, 1};
  elec_.Lattice.reset();

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
  elec_.createSK(); // needed by AoS J2 for ChiesaKEcorrection

  ParticleSetPool ptcl{c};
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  // make a ParticleSet Clone
  ParticleSet elec_clone(elec_);

  //diamondC_1x1x1
  const char* spo_xml = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_2x1x1.pwscf.h5\" tilematrix=\"2 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"2\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml);
  REQUIRE(okay);

  xmlNodePtr spo_root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(spo_root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  SPOSet* spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != nullptr);

  auto* det_up = new DiracDet(spo);
  det_up->set(0, 2);
  auto* det_dn = new DiracDet(spo);
  det_dn->set(2, 2);

  auto* slater_det = new SlaterDet(elec_);
  slater_det->add(det_up, 0);
  slater_det->add(det_dn, 1);

  TrialWaveFunction psi{c};
  psi.addComponent(slater_det, "SingleDet");

  const char* jas_input = "<tmp> \
<jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\"> \
   <correlation size=\"10\" speciesA=\"u\" speciesB=\"u\"> \
      <coefficients id=\"uu\" type=\"Array\"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201 -0.3253286875 -0.3624525145 -0.3958223107 -0.4268582166 -0.4394531176</coefficients> \
   </correlation> \
</jastrow> \
</tmp> \
";
  Libxml2Document doc_jas;
  okay = doc.parseFromString(jas_input);
  REQUIRE(okay);

  xmlNodePtr jas_root = doc.getRoot();
  xmlNodePtr jas1 = xmlFirstElementChild(jas_root);

  RadialJastrowBuilder jb(c, elec_);
  psi.addComponent(jb.buildComponent(jas1), "RadialJastrow");

#if !defined(QMC_CUDA)
  // initialize distance tables.
  elec_.update();
  double logpsi = psi.evaluateLog(elec_);

  //std::cout << "debug before YYY " << std::setprecision(16) << psi.getLogPsi() << " " << psi.getPhase()<< std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi == Approx(-4.546400274690973));
#else
  REQUIRE(logpsi == Approx(-5.932687817638513));
#endif

  // make a TrialWaveFunction Clone
  TrialWaveFunction* psi_clone = psi.makeClone(elec_clone);

  elec_clone.update();
  double logpsi_clone = psi_clone->evaluateLog(elec_clone);
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi_clone == Approx(-4.546400274690973));
#else
  REQUIRE(logpsi_clone == Approx(-5.932687817638513));
#endif

  const int moved_elec_id = 0;

  using PosType = QMCTraits::PosType;
  using RealType = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  PosType delta(0.1, 0.1, 0.2);

  elec_.makeMove(moved_elec_id, delta);

  ValueType r_all_val = psi.calcRatio(elec_, moved_elec_id);
  ValueType r_fermionic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::FERMIONIC);
  ValueType r_bosonic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::NONFERMIONIC);

  //std::cout << "debug YYY " << std::setprecision(16) << r_all_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_fermionic_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_bosonic_val << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(r_all_val == ComplexApprox(ValueType(0.1249724458593248, -1.301480654871811e-06)));
  REQUIRE(r_fermionic_val == ComplexApprox(ValueType(0.1363257116240675, -1.419715163773587e-06)));
#else
  REQUIRE(r_all_val == Approx(0.1249737473276797));
  REQUIRE(r_fermionic_val == ValueApprox(0.1363271313258139));
#endif
  REQUIRE(r_bosonic_val == ValueApprox(0.9167195562048454));

  psi.acceptMove(elec_, moved_elec_id);
  elec_.acceptMove(moved_elec_id);
  //std::cout << "debug before YYY getLogPsi " << std::setprecision(16) << psi.getLogPsi() << " " << psi.getPhase()<< std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(-6.626062273740937));
#else
  REQUIRE(psi.getLogPsi() == Approx(-8.012339402754488));
#endif

  elec_.update(true);
  psi.evaluateLog(elec_);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(-6.626062273740937));
#else
  REQUIRE(psi.getLogPsi() == Approx(-8.012339402754488));
#endif

  // testing batched interfaces
  std::vector<ParticleSet*> P_list(2, nullptr);
  P_list[0] = &elec_;
  P_list[1] = &elec_clone;

  std::vector<TrialWaveFunction*> WF_list(2, nullptr);
  WF_list[0] = &psi;
  WF_list[1] = psi_clone;

    //Temporary as switch to std::reference_wrapper proceeds
// testing batched interfaces
  RefVector<ParticleSet> p_ref_list{elec_,elec_clone};
  RefVector<TrialWaveFunction> wf_ref_list{psi, *psi_clone};

  elec_.flex_update(p_ref_list);
  psi.flex_evaluateLog(wf_ref_list, p_ref_list);
  //std::cout << "debug before YYY getLogPsi getPhase " << std::setprecision(16) << WF_list[0]->getLogPsi() << " " << WF_list[0]->getPhase()<< std::endl;
  //std::cout << "debug before YYY getLogPsi getPhase " << std::setprecision(16) << WF_list[1]->getLogPsi() << " " << WF_list[1]->getPhase()<< std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-6.626062273740933, -3.14160988592975)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-4.546400274690973, -3.141599471788856)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-8.012339402754462, 6.283185307179586)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-5.932687817638513, 6.283185307179586)));
#endif


  std::vector<GradType> grad_old(2);

  grad_old[0] = WF_list[0]->evalGrad(*P_list[0], moved_elec_id);
  grad_old[1] = WF_list[1]->evalGrad(*P_list[1], moved_elec_id);

  std::cout << "evalGrad " << std::setprecision(14)
            << grad_old[0][0] << " " << grad_old[0][1] << " " << grad_old[0][2] << " "
            << grad_old[1][0] << " " << grad_old[1][1] << " " << grad_old[1][2]
            << std::endl;  
  
  psi.flex_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);

#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(712.97995328773, 0.129153192215)));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(712.9798197679, 0.12906358823709)));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(-767.67252717682, -0.12986626645764)));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(117.94653318921, -0.011021468167199)));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(117.94652437244, -0.011021468209138)));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(-118.38325389331, 0.011021468183434)));
#else
  REQUIRE(grad_old[0][0] == Approx(712.8508040062));
  REQUIRE(grad_old[0][1] == Approx(712.85076008764));
  REQUIRE(grad_old[0][2] == Approx(-767.54266484264));
  REQUIRE(grad_old[1][0] == Approx(117.95755455321));
  REQUIRE(grad_old[1][1] == Approx(117.95754573648));
  REQUIRE(grad_old[1][2] == Approx(-118.39427525732));
#endif

  PosType delta_sign_changed(0.1, 0.1, -0.2);
  p_ref_list[0].get().makeMove(moved_elec_id, delta_sign_changed);
  p_ref_list[1].get().makeMove(moved_elec_id, delta_sign_changed);

  ValueType r_0 = wf_ref_list[0].get().calcRatio(p_ref_list[0].get(), moved_elec_id);
  GradType grad_temp;
  ValueType r_1 = wf_ref_list[1].get().calcRatioGrad(p_ref_list[1].get(), moved_elec_id, grad_temp);
  std::cout << "calcRatio calcRatioGrad " << std::setprecision(14)
            << r_0 << " " << r_1 << " "
            << grad_temp[0] << " " << grad_temp[1] << " " << grad_temp[2]
            << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(r_0 == ComplexApprox(ValueType(253.47871118162,-0.0016996449989914)));
  REQUIRE(r_1 == ComplexApprox(ValueType(36.908418548171,-0.001734416845066)));
  REQUIRE(grad_temp[0] == ComplexApprox(ValueType(1.4581181602054,0.00043036738827721)));
  REQUIRE(grad_temp[1] == ComplexApprox(ValueType(1.4581187980287,0.00043120067711778)));
  REQUIRE(grad_temp[2] == ComplexApprox(ValueType(-1.2913531933125,-0.00033355989998835)));
#else
  REQUIRE(r_0 == Approx(253.48041077515));
  REQUIRE(r_1 == Approx(36.910152948623));
  REQUIRE(grad_temp[0] == Approx(1.4576878373305));
  REQUIRE(grad_temp[1] == Approx(1.4576876419511));
  REQUIRE(grad_temp[2] == Approx(-1.291019667913));
#endif

  PosType delta_zero(0, 0, 0);
  p_ref_list[0].get().makeMove(moved_elec_id, delta_zero);
  p_ref_list[1].get().makeMove(moved_elec_id, delta);

  std::vector<PsiValueType> ratios(2);
  psi.flex_calcRatio(wf_ref_list, p_ref_list, moved_elec_id, ratios);
  std::cout << "calcRatio " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl;

#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(PsiValueType(1, 0)));
  REQUIRE(ratios[1] == ComplexApprox(PsiValueType(0.12497244585932,-1.3014806548718e-06)));
#else
  REQUIRE(ratios[0] == Approx(1));
  REQUIRE(ratios[1] == Approx(0.12497374732768));
#endif

  std::fill(ratios.begin(), ratios.end(), 0);
  std::vector<GradType> grad_new(2);

  ratios[0] = WF_list[0]->calcRatioGrad(*P_list[0], moved_elec_id, grad_new[0]);
  ratios[1] = WF_list[1]->calcRatioGrad(*P_list[1], moved_elec_id, grad_new[1]);

  std::cout << "calcRatioGrad " << std::setprecision(14)
            << ratios[0] << " " << ratios[1] << std::endl
            << grad_new[0][0] << " " << grad_new[0][1] << " " << grad_new[0][2] << " "
            << grad_new[1][0] << " " << grad_new[1][1] << " " << grad_new[1][2]
            << std::endl;

  //Temporary as switch to std::reference_wrapper proceeds
  // testing batched interfaces

  psi.flex_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id, ratios, grad_new);
#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(ValueType(1, 0)));
  REQUIRE(grad_new[0][0] == ComplexApprox(ValueType(712.97995328773,0.12915319218279)));
  REQUIRE(grad_new[0][1] == ComplexApprox(ValueType(712.9798197679,0.12906358820488)));
  REQUIRE(grad_new[0][2] == ComplexApprox(ValueType(-767.67252717682,-0.12986626642298)));
  REQUIRE(ratios[1] == ComplexApprox(ValueType(0.12497244585932,-1.3014806548718e-06)));
  REQUIRE(grad_new[1][0] == ComplexApprox(ValueType(712.97995328773,0.12915319219098)));
  REQUIRE(grad_new[1][1] == ComplexApprox(ValueType(712.9798197679,0.12906358821309)));
  REQUIRE(grad_new[1][2] == ComplexApprox(ValueType(-767.67252717682,-0.12986626643177)));
#else
  REQUIRE(ratios[0] == Approx(1));
  REQUIRE(grad_new[0][0] == Approx(712.85080400624));
  REQUIRE(grad_new[0][1] == Approx(712.85076008768));
  REQUIRE(grad_new[0][2] == Approx(-767.54266484268));
  REQUIRE(ratios[1] == Approx(0.12497374732768));
  REQUIRE(grad_new[1][0] == Approx(712.85080400622));
  REQUIRE(grad_new[1][1] == Approx(712.85076008766));
  REQUIRE(grad_new[1][2] == Approx(-767.54266484266));
#endif

  psi.flex_acceptMove(wf_ref_list, p_ref_list, moved_elec_id);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-6.626062273740933, -3.14160988592975)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-6.626062273740933, -3.14160988592975)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-8.012339402754462, 6.283185307179586)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-8.012339402754462, 6.283185307179586)));
#endif

  psi.flex_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);
#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(712.97995328773, 0.129153192215)));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(712.9798197679, 0.12906358823709)));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(-767.67252717682, -0.12986626645764)));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(712.97995328773, 0.129153192215)));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(712.9798197679, 0.12906358823709)));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(-767.67252717682, -0.12986626645764)));
#else
  REQUIRE(grad_old[0][0] == Approx(712.8508040062));
  REQUIRE(grad_old[0][1] == Approx(712.85076008764));
  REQUIRE(grad_old[0][2] == Approx(-767.54266484264));
  REQUIRE(grad_old[1][0] == Approx(712.8508040062));
  REQUIRE(grad_old[1][1] == Approx(712.85076008764));
  REQUIRE(grad_old[1][2] == Approx(-767.54266484264));
#endif


  //FIXME more thinking and fix about ownership and schope are needed for exiting clean
  delete psi_clone;
#endif

}

} // namespace qmcplusplus
