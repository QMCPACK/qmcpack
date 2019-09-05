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

using LogValueType = WaveFunctionComponent::LogValueType;

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

  auto* det_up = new DiracDeterminant<>(spo);
  det_up->set(0, 2);
  auto* det_dn = new DiracDeterminant<>(spo);
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

  RadialJastrowBuilder jb(elec_, psi);
  bool build_okay = jb.put(jas1);
  REQUIRE(build_okay);

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

  std::vector<ParticleSet*> P_list(2, nullptr);
  P_list[0] = &elec_;
  P_list[1] = &elec_clone;

  std::vector<TrialWaveFunction*> WF_list(2, nullptr);
  WF_list[0] = &psi;
  WF_list[1] = psi_clone;

  psi.flex_evaluateLog(WF_list, P_list);
#if defined(QMC_COMPLEX)
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == ComplexApprox(std::complex<RealType>(0.4351202455204972, 6.665972664860828)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == ComplexApprox(std::complex<RealType>(-0.1201465271523596, 6.345732826640545)));
#else
  REQUIRE(std::complex<RealType>(WF_list[0]->getLogPsi(), WF_list[0]->getPhase()) == ComplexApprox(std::complex<RealType>(-0.6365029797784554, 3.141592653589793)));
  REQUIRE(std::complex<RealType>(WF_list[1]->getLogPsi(), WF_list[1]->getPhase()) == ComplexApprox(std::complex<RealType>(-1.471840358291562, 3.141592653589793)));
#endif

  //FIXME more thinking and fix about ownership and schope are needed for exiting clean
  delete psi_clone;
#endif

}

} // namespace qmcplusplus
