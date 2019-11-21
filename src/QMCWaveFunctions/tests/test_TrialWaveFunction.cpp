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
#include "Particle/ParticleSet.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"

namespace qmcplusplus
{

#ifdef ENABLE_CUDA
using DiracDet = DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#else
using DiracDet = DiracDeterminant<DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

using LogValueType = TrialWaveFunction::LogValueType;
using PsiValueType = TrialWaveFunction::PsiValueType;
using GradType = TrialWaveFunction::GradType;

TEST_CASE("TrialWaveFunction", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion");
  ions_.create(2);
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  ions_.update();


  elec_.setName("elec");
  std::vector<int> ud(2);
  ud[0] = ud[1] = 2;
  elec_.create(ud);
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 1.0};
  elec_.R[2] = {1.0, 1.0, 0.0};
  elec_.R[3] = {1.0, 0.0, 1.0};

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
<determinantset type=\"einspline\" href=\"pwscf.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"2\"/> \
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
  REQUIRE(logpsi == Approx(-0.1201465271523596));
#else
  REQUIRE(logpsi == Approx(-1.471840358291562));
#endif

  // make a TrialWaveFunction Clone
  TrialWaveFunction* psi_clone = psi.makeClone(elec_clone);

  elec_clone.update();
  double logpsi_clone = psi_clone->evaluateLog(elec_clone);
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi_clone == Approx(-0.1201465271523596));
#else
  REQUIRE(logpsi_clone == Approx(-1.471840358291562));
#endif

  const int moved_elec_id = 0;

  using PosType = QMCTraits::PosType;
  using RealType = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  PosType delta(0.1, 0.1, 0.2);

  elec_.setActive(moved_elec_id);
  elec_.makeMove(moved_elec_id, delta);

  RealType r_val = psi.ratio(elec_, moved_elec_id);
  ValueType r_all_val = psi.calcRatio(elec_, moved_elec_id);
  ValueType r_fermionic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::FERMIONIC);
  ValueType r_bosonic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::NONFERMIONIC);

  //std::cout << "debug YYY " << std::setprecision(16) << r_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_all_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_fermionic_val << std::endl;
  //std::cout << "debug YYY " << std::setprecision(16) << r_bosonic_val << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(r_val == Approx(1.742405749016986));
  REQUIRE(r_all_val == ComplexApprox(ValueType(1.653821746120792, 0.5484992491019633)));
  REQUIRE(r_fermionic_val == ComplexApprox(ValueType(1.804065087219802, 0.598328295048828)));
#else
  REQUIRE(r_val == Approx(2.305591774210242));
  REQUIRE(r_all_val == Approx(2.305591774210242));
  REQUIRE(r_fermionic_val == ValueApprox(2.515045914101833));
#endif
  REQUIRE(r_bosonic_val == ValueApprox(0.9167195562048454));

  psi.acceptMove(elec_, moved_elec_id);
  elec_.acceptMove(moved_elec_id);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(0.4351202455204972));
#else
  REQUIRE(psi.getLogPsi() == Approx(-0.63650297977845492));
#endif

  elec_.update(true);
  psi.evaluateLog(elec_);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(0.4351202455204972));
#else
  REQUIRE(psi.getLogPsi() == Approx(-0.63650297977845492));
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
  
  psi.flex_evaluateLog(wf_ref_list, p_ref_list);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-0.1201465271523596, 6.345732826640545)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-1.471840358291562, 3.141592653589793)));
#endif

  P_list[0]->setActive(moved_elec_id);
  P_list[1]->setActive(moved_elec_id);

  std::vector<GradType> grad_old(2);

  grad_old[0] = WF_list[0]->evalGrad(*P_list[0], moved_elec_id);
  grad_old[1] = WF_list[1]->evalGrad(*P_list[1], moved_elec_id);

  std::cout << "evalGrad " << std::setprecision(14)
            << grad_old[0][0] << " " << grad_old[0][1] << " " << grad_old[0][2] << " "
            << grad_old[1][0] << " " << grad_old[1][1] << " " << grad_old[1][2]
            << std::endl;  
  
  psi.flex_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);
#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(47.387717528888, -8.7703065253151e-06)));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(-54.671696901113, -7.3126138879524)));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(6.6288917088321,7.3126230586018)));
#else
  REQUIRE(grad_old[0][0] == Approx(14.77249702264));
  REQUIRE(grad_old[0][1] == Approx(-20.385235323777));
  REQUIRE(grad_old[0][2] == Approx(4.8529516184558));
  REQUIRE(grad_old[1][0] == Approx(47.38770710732));
  REQUIRE(grad_old[1][1] == Approx(-63.361119579044));
  REQUIRE(grad_old[1][2] == Approx(15.318325284049));
#endif

  PosType delta_zero(0, 0, 0);
  p_ref_list[0].get().makeMove(moved_elec_id, delta_zero);
  p_ref_list[1].get().makeMove(moved_elec_id, delta);

  std::vector<PsiValueType> ratios(2);
  psi.flex_calcRatio(wf_ref_list, p_ref_list, moved_elec_id, ratios);
  std::cout << "calcRatio " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(PsiValueType(1, 0)));
  REQUIRE(ratios[1] == ComplexApprox(PsiValueType(1.6538214581548,0.54849918598717)));
#else
  REQUIRE(ratios[0] == Approx(1));
  REQUIRE(ratios[1] == Approx(2.3055913093424));
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
  
  psi.flex_ratioGrad(wf_ref_list, p_ref_list, moved_elec_id, ratios, grad_new);
#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(ValueType(1, 0)));
  REQUIRE(grad_new[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  REQUIRE(grad_new[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  REQUIRE(grad_new[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  REQUIRE(ratios[1] == ComplexApprox(ValueType(1.6538214581548,0.54849918598717)));
  REQUIRE(grad_new[1][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  REQUIRE(grad_new[1][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  REQUIRE(grad_new[1][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
#else
  REQUIRE(ratios[0] == Approx(1));
  REQUIRE(grad_new[0][0] == Approx(14.77249702264));
  REQUIRE(grad_new[0][1] == Approx(-20.385235323777));
  REQUIRE(grad_new[0][2] == Approx(4.8529516184558));
  REQUIRE(ratios[1] == Approx(2.3055913093424));
  REQUIRE(grad_new[1][0] == Approx(14.77249702264));
  REQUIRE(grad_new[1][1] == Approx(-20.385235323777));
  REQUIRE(grad_new[1][2] == Approx(4.8529516184558));
#endif

  psi.flex_acceptMove(wf_ref_list, p_ref_list, moved_elec_id);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == LogComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
#endif

  psi.flex_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);
#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(18.817970466022, -6.5837500306076)));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(-22.840838391977, 3.9963373883645)));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(3.8805320617146, 1.5825508129169)));
#else
  REQUIRE(grad_old[0][0] == Approx(14.77249702264));
  REQUIRE(grad_old[0][1] == Approx(-20.385235323777));
  REQUIRE(grad_old[0][2] == Approx(4.8529516184558));
  REQUIRE(grad_old[1][0] == Approx(14.77249702264));
  REQUIRE(grad_old[1][1] == Approx(-20.385235323777));
  REQUIRE(grad_old[1][2] == Approx(4.8529516184558));
#endif


  //FIXME more thinking and fix about ownership and schope are needed for exiting clean
  delete psi_clone;
#endif

}

} // namespace qmcplusplus
