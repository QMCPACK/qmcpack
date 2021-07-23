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
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"

#include <regex>

namespace qmcplusplus
{
using LogValueType = TrialWaveFunction::LogValueType;
using PsiValueType = TrialWaveFunction::PsiValueType;
using GradType     = TrialWaveFunction::GradType;

struct double_tag
{};
struct float_tag
{};

template<class DiracDet, class SPO_precision>
void testTrialWaveFunction_diamondC_2x1x1(const int ndelay)
{
#if defined(MIXED_PRECISION)
  const double precision = 1e-4;
#else
  const double precision = std::is_same<SPO_precision, float_tag>::value ? 1e-4 : 1e-8;
#endif
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  auto ions_uptr = std::make_unique<ParticleSet>();
  auto elec_uptr = std::make_unique<ParticleSet>();
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

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
  elec_.Lattice.R(0, 0)   = 6.7463223;
  elec_.Lattice.R(0, 1)   = 6.7463223;
  elec_.Lattice.R(0, 2)   = 0.0;
  elec_.Lattice.R(1, 0)   = 0.0;
  elec_.Lattice.R(1, 1)   = 3.37316115;
  elec_.Lattice.R(1, 2)   = 3.37316115;
  elec_.Lattice.R(2, 0)   = 3.37316115;
  elec_.Lattice.R(2, 1)   = 0.0;
  elec_.Lattice.R(2, 2)   = 3.37316115;
  elec_.Lattice.BoxBConds = {1, 1, 1};
  elec_.Lattice.reset();

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec_.addTable(ions_);
  elec_.resetGroups();
  elec_.createSK(); // needed by AoS J2 for ChiesaKEcorrection

  ParticleSetPool ptcl{c};
  ptcl.addParticleSet(std::move(elec_uptr));
  ptcl.addParticleSet(std::move(ions_uptr));

  // make a ParticleSet Clone
  ParticleSet elec_clone(elec_);

  //diamondC_1x1x1
  std::string spo_xml = "<tmp> \
               <determinantset type=\"einspline\" href=\"diamondC_2x1x1.pwscf.h5\" tilematrix=\"2 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.5\" precision=\"float\" size=\"2\"/> \
               </tmp> \
               ";
#ifndef MIXED_PRECISION
  if (std::is_same<SPO_precision, double_tag>::value)
  {
    spo_xml = std::regex_replace(spo_xml, std::regex("float"), "double");
  }
#endif
  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml);
  REQUIRE(okay);

  xmlNodePtr spo_root = doc.getRoot();
  xmlNodePtr ein1     = xmlFirstElementChild(spo_root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != nullptr);

  auto* det_up = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
  det_up->set(0, 2, ndelay);
  auto* det_dn = new DiracDet(std::unique_ptr<SPOSet>(spo->makeClone()));
  det_dn->set(2, 2, ndelay);

  auto slater_det = std::make_unique<SlaterDet>(elec_);
  slater_det->add(det_up, 0);
  slater_det->add(det_dn, 1);

  TrialWaveFunction psi;
  psi.addComponent(std::move(slater_det));

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
  xmlNodePtr jas1     = xmlFirstElementChild(jas_root);

  RadialJastrowBuilder jb(c, elec_);
  psi.addComponent(jb.buildComponent(jas1));

  ResourceCollection res_col("test_determinant");
  psi.createResource(res_col);
  psi.acquireResource(res_col);

#if !defined(QMC_CUDA)
  // initialize distance tables.
  elec_.update();
  double logpsi = psi.evaluateLog(elec_);

  std::cout << "debug before YYY logpsi " << std::setprecision(16) << psi.getLogPsi() << " " << psi.getPhase()
            << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi == Approx(-4.546410485374186));
#else
  REQUIRE(logpsi == Approx(-5.932711221043984));
#endif

  auto logpsi_cplx = psi.evaluateGL(elec_, false);
#if defined(QMC_COMPLEX)
  REQUIRE(std::real(logpsi_cplx) == Approx(-4.546410485374186));
#else
  REQUIRE(std::real(logpsi_cplx) == Approx(-5.932711221043984));
#endif

  logpsi_cplx = psi.evaluateGL(elec_, true);
#if defined(QMC_COMPLEX)
  REQUIRE(std::real(logpsi_cplx) == Approx(-4.546410485374186));
#else
  REQUIRE(std::real(logpsi_cplx) == Approx(-5.932711221043984));
#endif

  // make a TrialWaveFunction Clone
  std::unique_ptr<TrialWaveFunction> psi_clone(psi.makeClone(elec_clone));

  elec_clone.update();

  res_col.rewind();
  psi.releaseResource(res_col);
  res_col.rewind();
  psi_clone->acquireResource(res_col);
  
  
  double logpsi_clone = psi_clone->evaluateLog(elec_clone);
#if defined(QMC_COMPLEX)
  REQUIRE(logpsi_clone == Approx(-4.546410485374186));
#else
  REQUIRE(logpsi_clone == Approx(-5.932711221043984));
#endif

  const int moved_elec_id = 0;

  using PosType   = QMCTraits::PosType;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  PosType delta(0.1, 0.1, 0.2);

  elec_.makeMove(moved_elec_id, delta);

  res_col.rewind();
  psi_clone->releaseResource(res_col);
  res_col.rewind();
  psi.acquireResource(res_col);

  ValueType r_all_val       = psi.calcRatio(elec_, moved_elec_id);
  ValueType r_fermionic_val = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::FERMIONIC);
  ValueType r_bosonic_val   = psi.calcRatio(elec_, moved_elec_id, TrialWaveFunction::ComputeType::NONFERMIONIC);

  std::cout << "YYY r_all_val " << std::setprecision(16) << r_all_val << std::endl;
  std::cout << "YYY r_fermionic_val " << std::setprecision(16) << r_fermionic_val << std::endl;
  std::cout << "YYY r_bosonic_val " << std::setprecision(16) << r_bosonic_val << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(r_all_val == ComplexApprox(std::complex<RealType>(0.1248738460467855, 0)).epsilon(2e-5));
  REQUIRE(r_fermionic_val == ComplexApprox(std::complex<RealType>(0.1362181543980086, 0)).epsilon(2e-5));
#else
  REQUIRE(r_all_val == Approx(0.1248738460469678));
  REQUIRE(r_fermionic_val == ValueApprox(0.1362181543982075));
#endif
  REQUIRE(r_bosonic_val == ValueApprox(0.9167195562048454));

  psi.acceptMove(elec_, moved_elec_id);
  elec_.acceptMove(moved_elec_id);
  std::cout << "before YYY getLogPsi " << std::setprecision(16) << psi.getLogPsi() << " " << psi.getPhase()
            << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(-6.626861768296886).epsilon(5e-5));
#else
  REQUIRE(psi.getLogPsi() == Approx(-8.013162503965223));
#endif

  elec_.update(true);
  psi.evaluateLog(elec_);
#if defined(QMC_COMPLEX)
  REQUIRE(psi.getLogPsi() == Approx(-6.626861768296886).epsilon(5e-5));
#else
  REQUIRE(psi.getLogPsi() == Approx(-8.013162503965223));
#endif

  // testing batched interfaces
  std::vector<ParticleSet*> P_list(2, nullptr);
  P_list[0] = &elec_;
  P_list[1] = &elec_clone;

  std::vector<TrialWaveFunction*> WF_list(2, nullptr);
  WF_list[0] = &psi;
  WF_list[1] = psi_clone.get();

  //Temporary as switch to std::reference_wrapper proceeds
  // testing batched interfaces
  RefVectorWithLeader<ParticleSet> p_ref_list(elec_, {elec_, elec_clone});
  RefVectorWithLeader<TrialWaveFunction> wf_ref_list(psi, {psi, *psi_clone});

  ParticleSet::mw_update(p_ref_list);
  TrialWaveFunction::mw_evaluateLog(wf_ref_list, p_ref_list);
  std::cout << "before YYY [0] getLogPsi getPhase " << std::setprecision(16) << WF_list[0]->getLogPsi() << " "
            << WF_list[0]->getPhase() << std::endl;
  std::cout << "before YYY [1] getLogPsi getPhase " << std::setprecision(16) << WF_list[1]->getLogPsi() << " "
            << WF_list[1]->getPhase() << std::endl;
#if defined(QMC_COMPLEX)
  CHECK(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-6.626861768296848, -3.141586279082042)));
  CHECK(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-4.546410485374186, -3.141586279080522)));
#else
  CHECK(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-8.013162503965042, 6.283185307179586)));
  CHECK(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-5.932711221043984, 6.283185307179586)));
#endif


  std::vector<GradType> grad_old(2);

  grad_old[0] = WF_list[0]->evalGrad(*P_list[0], moved_elec_id);
  grad_old[1] = WF_list[1]->evalGrad(*P_list[1], moved_elec_id);

  std::cout << "evalGrad " << std::setprecision(14) << grad_old[0][0] << " " << grad_old[0][1] << " " << grad_old[0][2]
            << " " << grad_old[1][0] << " " << grad_old[1][1] << " " << grad_old[1][2] << std::endl;

  TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id, grad_old);

#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(713.71203320653, 0.020838031926441)).epsilon(precision));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(713.71203320654, 0.020838031928415)).epsilon(precision));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(-768.42842826889, -0.020838032018344)).epsilon(precision));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(118.02653358655, -0.0022419843505538)).epsilon(precision));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(118.02653358655, -0.0022419843498631)).epsilon(precision));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(-118.46325895634, 0.0022419843493758)).epsilon(precision));
#else
  REQUIRE(grad_old[0][0] == Approx(713.69119517454).epsilon(2. * precision));
  REQUIRE(grad_old[0][1] == Approx(713.69119517455).epsilon(2. * precision));
  REQUIRE(grad_old[0][2] == Approx(-768.40759023681).epsilon(2. * precision));
  REQUIRE(grad_old[1][0] == Approx(118.0287755709).epsilon(precision));
  REQUIRE(grad_old[1][1] == Approx(118.0287755709).epsilon(precision));
  REQUIRE(grad_old[1][2] == Approx(-118.46550094069).epsilon(precision));
#endif
  PosType delta_sign_changed(0.1, 0.1, -0.2);
  p_ref_list[0].makeMove(moved_elec_id, delta_sign_changed);
  p_ref_list[1].makeMove(moved_elec_id, delta_sign_changed);

  ValueType r_0 = wf_ref_list[0].calcRatio(p_ref_list[0], moved_elec_id);
  GradType grad_temp;
  ValueType r_1 = wf_ref_list[1].calcRatioGrad(p_ref_list[1], moved_elec_id, grad_temp);
  std::cout << "calcRatio calcRatioGrad " << std::setprecision(14) << r_0 << " " << r_1 << " " << grad_temp[0] << " "
            << grad_temp[1] << " " << grad_temp[2] << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(r_0 == ComplexApprox(ValueType(253.71869245791, -0.00034808849808193)).epsilon(1e-4));
  REQUIRE(r_1 == ComplexApprox(ValueType(36.915636007059, -6.4240180082292e-05)).epsilon(1e-5));
  REQUIRE(grad_temp[0] == ComplexApprox(ValueType(1.4567170375539, 0.00027263382943948)));
  REQUIRE(grad_temp[1] == ComplexApprox(ValueType(1.4567170375539, 0.00027263382945093)));
  REQUIRE(grad_temp[2] == ComplexApprox(ValueType(-1.2930978490431, -0.00027378452214318)));
#else
  REQUIRE(r_0 == Approx(253.71904054638).epsilon(2e-4));
  REQUIRE(r_1 == Approx(36.915700247239));
  REQUIRE(grad_temp[0] == Approx(1.4564444046733));
  REQUIRE(grad_temp[1] == Approx(1.4564444046734));
  REQUIRE(grad_temp[2] == Approx(-1.2928240654738));
#endif

  PosType delta_zero(0, 0, 0);
  p_ref_list[0].makeMove(moved_elec_id, delta_zero);
  p_ref_list[1].makeMove(moved_elec_id, delta);

  std::vector<PsiValueType> ratios(2);
  TrialWaveFunction::mw_calcRatio(wf_ref_list, p_ref_list, moved_elec_id, ratios);
  std::cout << "mixed move calcRatio " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl;

#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(PsiValueType(1, 0)).epsilon(5e-5));
#if defined(MIXED_PRECISION)
  REQUIRE(ratios[1] == ComplexApprox(PsiValueType(0.12487384604679, 0)).epsilon(2e-5));
#else
  REQUIRE(ratios[1] == ComplexApprox(PsiValueType(0.12487384604679, 0)));
#endif
#else
  REQUIRE(ratios[0] == Approx(1).epsilon(5e-5));
  REQUIRE(ratios[1] == Approx(0.12487384604697));
#endif

  std::fill(ratios.begin(), ratios.end(), 0);
  std::vector<GradType> grad_new(2);

  ratios[0] = WF_list[0]->calcRatioGrad(*P_list[0], moved_elec_id, grad_new[0]);
  ratios[1] = WF_list[1]->calcRatioGrad(*P_list[1], moved_elec_id, grad_new[1]);

  std::cout << "calcRatioGrad " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl
            << grad_new[0][0] << " " << grad_new[0][1] << " " << grad_new[0][2] << " " << grad_new[1][0] << " "
            << grad_new[1][1] << " " << grad_new[1][2] << std::endl;

  //Temporary as switch to std::reference_wrapper proceeds
  // testing batched interfaces

  TrialWaveFunction::mw_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id, ratios, grad_new);
  std::cout << "flex_calcRatioGrad " << std::setprecision(14) << ratios[0] << " " << ratios[1] << std::endl
            << grad_new[0][0] << " " << grad_new[0][1] << " " << grad_new[0][2] << " " << grad_new[1][0] << " "
            << grad_new[1][1] << " " << grad_new[1][2] << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(ratios[0] == ComplexApprox(ValueType(1, 0)).epsilon(5e-5));
  REQUIRE(grad_new[0][0] == ComplexApprox(ValueType(713.71203320653, 0.020838031942702)).epsilon(precision));
  REQUIRE(grad_new[0][1] == ComplexApprox(ValueType(713.71203320654, 0.020838031944677)).epsilon(precision));
  REQUIRE(grad_new[0][2] == ComplexApprox(ValueType(-768.42842826889, -0.020838032035842)).epsilon(precision));
  REQUIRE(ratios[1] == ComplexApprox(ValueType(0.12487384604679, 0)));
  REQUIRE(grad_new[1][0] == ComplexApprox(ValueType(713.71203320656, 0.020838031892613)).epsilon(precision));
  REQUIRE(grad_new[1][1] == ComplexApprox(ValueType(713.71203320657, 0.020838031894628)).epsilon(precision));
  REQUIRE(grad_new[1][2] == ComplexApprox(ValueType(-768.42842826892, -0.020838031981896)).epsilon(precision));
#else
  REQUIRE(ratios[0] == Approx(1).epsilon(5e-5));
  REQUIRE(grad_new[0][0] == Approx(713.69119517463).epsilon(precision));
  REQUIRE(grad_new[0][1] == Approx(713.69119517463).epsilon(precision));
  REQUIRE(grad_new[0][2] == Approx(-768.40759023689).epsilon(precision));
  REQUIRE(ratios[1] == Approx(0.12487384604697));
  REQUIRE(grad_new[1][0] == Approx(713.69119517467).epsilon(precision));
  REQUIRE(grad_new[1][1] == Approx(713.69119517468).epsilon(precision));
  REQUIRE(grad_new[1][2] == Approx(-768.40759023695).epsilon(precision));
#endif

  std::vector<bool> isAccepted(2, true);
  TrialWaveFunction::mw_accept_rejectMove(wf_ref_list, p_ref_list, moved_elec_id, isAccepted, true);
  std::cout << "flex_acceptMove WF_list[0] getLogPsi getPhase " << std::setprecision(16) << WF_list[0]->getLogPsi()
            << " " << WF_list[0]->getPhase() << std::endl;
  std::cout << "flex_acceptMove WF_list[1] getLogPsi getPhase " << std::setprecision(16) << WF_list[1]->getLogPsi()
            << " " << WF_list[1]->getPhase() << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-6.626861768296848, -3.141586279082065)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-6.626861768296886, -3.141586279081995)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-8.013162503965155, 6.283185307179586)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) ==
          LogComplexApprox(std::complex<RealType>(-8.013162503965223, 6.283185307179586)));
#endif

  ParticleSet::mw_accept_rejectMove(p_ref_list, moved_elec_id, isAccepted, false);

  const int moved_elec_id_next = 1;
  TrialWaveFunction::mw_evalGrad(wf_ref_list, p_ref_list, moved_elec_id_next, grad_old);
  std::cout << "evalGrad next electron " << std::setprecision(14) << grad_old[0][0] << " " << grad_old[0][1] << " "
            << grad_old[0][2] << " " << grad_old[1][0] << " " << grad_old[1][1] << " " << grad_old[1][2] << std::endl;
#if defined(MIXED_PRECISION)
#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(-114.82740072726, -7.605305979232e-05)).epsilon(precision));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(-93.980772428401, -7.605302517238e-05)).epsilon(precision));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(64.050803536571, 7.6052975324197e-05)).epsilon(precision));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(-114.82740072726, -7.605305979232e-05)).epsilon(precision));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(-93.980772428401, -7.605302517238e-05)).epsilon(precision));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(64.050803536571, 7.6052975324197e-05)).epsilon(precision));
#else
  REQUIRE(grad_old[0][0] == Approx(-114.82732467419).epsilon(precision));
  REQUIRE(grad_old[0][1] == Approx(-93.98069637537).epsilon(precision));
  REQUIRE(grad_old[0][2] == Approx(64.050727483593).epsilon(precision));
  REQUIRE(grad_old[1][0] == Approx(-114.82732467419).epsilon(precision));
  REQUIRE(grad_old[1][1] == Approx(-93.98069637537).epsilon(precision));
  REQUIRE(grad_old[1][2] == Approx(64.050727483593).epsilon(precision));
#endif
#else
#if defined(QMC_COMPLEX)
  REQUIRE(grad_old[0][0] == ComplexApprox(ValueType(-114.82740072726, -7.605305979232e-05)).epsilon(precision));
  REQUIRE(grad_old[0][1] == ComplexApprox(ValueType(-93.980772428401, -7.605302517238e-05)).epsilon(precision));
  REQUIRE(grad_old[0][2] == ComplexApprox(ValueType(64.050803536571, 7.6052975324197e-05)).epsilon(precision));
  REQUIRE(grad_old[1][0] == ComplexApprox(ValueType(-114.82740072726, -7.605305979232e-05)).epsilon(precision));
  REQUIRE(grad_old[1][1] == ComplexApprox(ValueType(-93.980772428401, -7.605302517238e-05)).epsilon(precision));
  REQUIRE(grad_old[1][2] == ComplexApprox(ValueType(64.050803536571, 7.6052975324197e-05)).epsilon(precision));
#else
  REQUIRE(grad_old[0][0] == Approx(-114.82732467419).epsilon(precision));
  REQUIRE(grad_old[0][1] == Approx(-93.98069637537).epsilon(precision));
  REQUIRE(grad_old[0][2] == Approx(64.050727483593).epsilon(precision));
  REQUIRE(grad_old[1][0] == Approx(-114.82732467419).epsilon(precision));
  REQUIRE(grad_old[1][1] == Approx(-93.98069637537).epsilon(precision));
  REQUIRE(grad_old[1][2] == Approx(64.050727483593).epsilon(precision));
#endif
#endif

  std::vector<PosType> displ(2);
  displ[0] = displ[1] = {0.1, 0.2, 0.3};

  ParticleSet::mw_makeMove(p_ref_list, moved_elec_id_next, displ);
  TrialWaveFunction::mw_calcRatioGrad(wf_ref_list, p_ref_list, moved_elec_id_next, ratios, grad_new);
  std::cout << "ratioGrad next electron " << std::setprecision(14) << grad_new[0][0] << " " << grad_new[0][1] << " "
            << grad_new[0][2] << " " << grad_new[1][0] << " " << grad_new[1][1] << " " << grad_new[1][2] << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(grad_new[0][0] == ComplexApprox(ValueType(9.6073058494562, -1.4375146770852e-05)).epsilon(8e-5));
  REQUIRE(grad_new[0][1] == ComplexApprox(ValueType(6.3111018321898, -1.4375146510386e-05)).epsilon(8e-5));
  REQUIRE(grad_new[0][2] == ComplexApprox(ValueType(-3.2027658046121, 1.4375146020225e-05)).epsilon(8e-5));
  REQUIRE(grad_new[1][0] == ComplexApprox(ValueType(9.6073058494562, -1.4375146770852e-05)).epsilon(8e-5));
  REQUIRE(grad_new[1][1] == ComplexApprox(ValueType(6.3111018321898, -1.4375146510386e-05)).epsilon(8e-5));
  REQUIRE(grad_new[1][2] == ComplexApprox(ValueType(-3.2027658046121, 1.4375146020225e-05)).epsilon(8e-5));
#else
  REQUIRE(grad_new[0][0] == Approx(9.607320224603).epsilon(1e-4));
  REQUIRE(grad_new[0][1] == Approx(6.3111162073363).epsilon(1e-4));
  REQUIRE(grad_new[0][2] == Approx(-3.2027801797581).epsilon(1e-4));
  REQUIRE(grad_new[1][0] == Approx(9.607320224603).epsilon(1e-4));
  REQUIRE(grad_new[1][1] == Approx(6.3111162073363).epsilon(1e-4));
  REQUIRE(grad_new[1][2] == Approx(-3.2027801797581).epsilon(1e-4));
#endif

  isAccepted[0] = true;
  isAccepted[1] = false;
  TrialWaveFunction::mw_accept_rejectMove(wf_ref_list, p_ref_list, moved_elec_id_next, isAccepted, true);

  ParticleSet::mw_accept_rejectMove(p_ref_list, moved_elec_id_next, isAccepted, false);
  TrialWaveFunction::mw_completeUpdates(wf_ref_list);
  TrialWaveFunction::mw_evaluateGL(wf_ref_list, p_ref_list, false);

  
  std::cout << "invMat next electron " << std::setprecision(14) << det_up->get_det_engine().get_psiMinv()[0][0] << " "
            << det_up->get_det_engine().get_psiMinv()[0][1] << " " << det_up->get_det_engine().get_psiMinv()[1][0] << " " << det_up->get_det_engine().get_psiMinv()[1][1]
            << " " << std::endl;
#if defined(QMC_COMPLEX)
  REQUIRE(det_up->get_det_engine().get_psiMinv()[0][0] == ComplexApprox(ValueType(38.503358805635, -38.503358805645)).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[0][1] == ComplexApprox(ValueType(-31.465077529568, 31.465077529576)).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[1][0] == ComplexApprox(ValueType(-27.188228530061, 27.188228530068)).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[1][1] == ComplexApprox(ValueType(22.759962501254, -22.75996250126)).epsilon(1e-4));
#else
  REQUIRE(det_up->get_det_engine().get_psiMinv()[0][0] == Approx(77.0067176113).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[0][1] == Approx(-62.9301550592).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[1][0] == Approx(-54.376457060136).epsilon(1e-4));
  REQUIRE(det_up->get_det_engine().get_psiMinv()[1][1] == Approx(45.51992500251).epsilon(1e-4));
#endif
  std::vector<LogValueType> log_values(wf_ref_list.size());
  TrialWaveFunction::mw_evaluateGL(wf_ref_list, p_ref_list, false);
  for (int iw = 0; iw < log_values.size(); iw++)
    log_values[iw] = {wf_ref_list[iw].getLogPsi(), wf_ref_list[iw].getPhase()};
  TrialWaveFunction::mw_evaluateGL(wf_ref_list, p_ref_list, true);
  for (int iw = 0; iw < log_values.size(); iw++)
    REQUIRE(LogComplexApprox(log_values[iw]) == LogValueType{wf_ref_list[iw].getLogPsi(), wf_ref_list[iw].getPhase()});

#endif
}

TEST_CASE("TrialWaveFunction_diamondC_2x1x1", "[wavefunction]")
{
  using VT   = QMCTraits::ValueType;
  using FPVT = QMCTraits::QTFull::ValueType;


#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>, float_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>, float_tag>(2);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>, double_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixDelayedUpdateCUDA<VT, FPVT>>, double_tag>(2);
#endif
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixUpdateOMPTarget<VT, FPVT>>, float_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixUpdateOMPTarget<VT, FPVT>>, float_tag>(2);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixUpdateOMPTarget<VT, FPVT>>, double_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminantBatched<MatrixUpdateOMPTarget<VT, FPVT>>, double_tag>(2);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminant<DelayedUpdate<VT, FPVT>>, float_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminant<DelayedUpdate<VT, FPVT>>, float_tag>(2);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminant<DelayedUpdate<VT, FPVT>>, double_tag>(1);
  testTrialWaveFunction_diamondC_2x1x1<DiracDeterminant<DelayedUpdate<VT, FPVT>>, double_tag>(2);
}

} // namespace qmcplusplus
