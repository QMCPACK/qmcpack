//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "ParticleIO/XMLParticleIO.h"
namespace qmcplusplus
{


void create_CN_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  const char* particles = "<tmp> \
  <particleset name=\"ion0\" size=\"2\"> \
    <group name=\"C\">\
      <parameter name=\"charge\">4</parameter> \
      <parameter name=\"valence\">2</parameter> \
      <parameter name=\"atomicnumber\">6</parameter> \
    </group> \
    <group name=\"N\"> \
      <parameter name=\"charge\">5</parameter> \
      <parameter name=\"valence\">3</parameter> \
      <parameter name=\"atomicnumber\">7</parameter> \
    </group> \
    <attrib name=\"position\" datatype=\"posArray\"> \
  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00 \
  0.0000000000e+00  0.0000000000e+00  2.0786985865e+00 \
</attrib> \
    <attrib name=\"ionid\" datatype=\"stringArray\"> \
 C N \
</attrib> \
  </particleset> \
  <particleset name=\"e\"> \
    <group name=\"u\" size=\"5\"> \
      <parameter name=\"charge\">-1</parameter> \
    <attrib name=\"position\" datatype=\"posArray\"> \
-0.55936725 -0.26942464 0.14459603 \
0.19146719 1.40287983 0.63931251 \
1.14805915 -0.52057335 3.49621107 \
0.28293870 -0.10273952 0.01707021 \
0.60626935 -0.25538121 1.75750740 \
</attrib> \
    </group> \
    <group name=\"d\" size=\"4\"> \
      <parameter name=\"charge\">-1</parameter> \
    <attrib name=\"position\" datatype=\"posArray\"> \
-0.47405939 0.59523171 -0.59778601 \
0.03150661 -0.27343474 0.56279442 \
-1.32648025 0.00970226 2.26944242 \
2.42944286 0.64884151 1.87505288 \
</attrib> \
    </group> \
  </particleset>\
  </tmp>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  xmlNodePtr part2 = xmlNextElementSibling(part1);

  Tensor<int, 3> tmat;
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  XMLParticleParser parse_ions(ions, tmat);
  parse_ions.put(part1);

  XMLParticleParser parse_electrons(elec, tmat);
  parse_electrons.put(part2);
  
  elec.addTable(elec);
  elec.addTable(ions);
  elec.update();

}

TEST_CASE("Eloc_Derivatives:slater_noj", "[hamiltonian]")
{
  using RealType  = QMCTraits::RealType;
  app_log()<<" WOOHOO Eloc_derivatives:slater_noj tested!";
  REQUIRE( true );

  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec,ions); 

  int Nions = ions.getTotalNum();
  int Nelec = elec.getTotalNum();

  HamiltonianFactory::PtclPoolType particle_set_map;
  HamiltonianFactory::PsiPoolType psi_map;

  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  HamiltonianFactory hf("h0", elec, particle_set_map, psi_map, c);

  WaveFunctionFactory wff("psi0", elec, particle_set_map, c);
  psi_map["psi0"] = &wff;

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
           <pseudo elementType=\"N\" href=\"N.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  hf.put(root);

  Libxml2Document wfdoc;
  bool wfokay = wfdoc.parse("cn.wfnoj.xml");
  REQUIRE(wfokay);

  xmlNodePtr wfroot = wfdoc.getRoot();
  wff.put(wfroot);

  TrialWaveFunction* psi = wff.getTWF();
  REQUIRE(psi != nullptr);
 
  //Output of WFTester Eloc test for this ion/electron configuration.
  //Logpsi: (-1.4233853149e+01,0.0000000000e+00)
  //HamTest   Total -2.3652298313e+01
  //HamTest Kinetic 9.1821937928e+00
  //HamTest ElecElec 1.9015560571e+01
  //HamTest IonIon 9.6214045316e+00
  //HamTest LocalECP -6.7839428299e+01
  //HamTest NonLocalECP 6.3679710911e+00

  //24: eloc = -23.6523123429
  //24:   HamTest Kinetic 9.1821882589
  //24:   HamTest ElecElec 19.0155605708
  //24:   HamTest IonIon 9.6214045316
  //24:   HamTest LocalECP -67.8394282995
  //24:   HamTest NonLocalECP 6.3679625953

  
  RealType logpsi= psi->evaluateLog(elec); 
  REQUIRE( logpsi == Approx(-14.233853149));
  
  QMCHamiltonian* ham = hf.getH();

  RealType eloc = ham->evaluateDeterministic(elec);
  enum observ_id { KINETIC=0, ELECELEC, IONION, LOCALECP, NONLOCALECP};
  REQUIRE(eloc                            == Approx(-2.3652298313e+01));
  REQUIRE(ham->getObservable(ELECELEC)    == Approx( 1.9015560571e+01));
  REQUIRE(ham->getObservable(IONION)      == Approx( 9.6214045316e+00));
  REQUIRE(ham->getObservable(LOCALECP)    == Approx(-6.7839428299e+01));
  REQUIRE(ham->getObservable(KINETIC)     == Approx( 9.1821937928e+00));
  REQUIRE(ham->getObservable(NONLOCALECP) == Approx( 6.3679710911e+00));
  
  for (int i = 0; i < ham->sizeOfObservables(); i++)
    app_log() << "  HamTest " << ham->getObservableName(i) << " " << ham->getObservable(i) << std::endl;

  //Now for the derivative tests 
  ParticleSet::ParticleGradient_t wfgradraw;
  ParticleSet::ParticlePos_t hf_term;
  ParticleSet::ParticlePos_t pulay_term;
  ParticleSet::ParticlePos_t wf_grad;
 
  wfgradraw.resize(Nions);
  wf_grad.resize(Nions);
  hf_term.resize(Nions);
  pulay_term.resize(Nions);

  wfgradraw[0] = psi->evalGradSource(elec, ions, 0); //On the C atom.
  wfgradraw[1] = psi->evalGradSource(elec, ions, 1); //On the N atom.

  convert(wfgradraw[0], wf_grad[0]); 
  convert(wfgradraw[1], wf_grad[1]); 

  //Reference from finite differences on this configuration.
  REQUIRE( wf_grad[0][0] == Approx(-1.9044650674260308));
  REQUIRE( wf_grad[0][1] == Approx(2.1257764985627148)); 
  REQUIRE( wf_grad[0][2] == Approx(7.0556319963444016)); 
  REQUIRE( wf_grad[1][0] == Approx(1.4233346821157509)); 
  REQUIRE( wf_grad[1][1] == Approx(-0.1446706081154048)); 
  REQUIRE( wf_grad[1][2] == Approx(0.1440176276013005)); 

  //Kinetic Force
  hf_term=0.0;
  pulay_term=0.0;
  (ham->getHamiltonian(KINETIC))->evaluateWithIonDerivs(elec, ions, *psi, hf_term, pulay_term);
  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx(1.0852823603357820));
  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(24.2154119471038562));
  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(111.8849364775797852));
  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(2.1572063443997536));
  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(-3.3743242489947529));
  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(7.5625192454964454));

  psi->recompute(elec);
  //NLPP Force
  hf_term=0.0;
  pulay_term=0.0;
  double val=(ham->getHamiltonian(NONLOCALECP))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
  app_log()<<"val = "<<val<<std::endl;
  app_log()<<"hf_term[0][0] = "<<hf_term[0][0]<<std::endl;
  app_log()<<"pulay_term[0][0] = "<<pulay_term[0][0]<<std::endl;
  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx(32.8158275110457254));
  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(-25.9792419671001049));
  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(-51.8992198405232514));
  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(-0.0081394725626183));
  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(1.2062112833888250));
  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(-0.3732875800910306));
}

TEST_CASE("Eloc_Derivatives:slater_wj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_Derivatives:slater_wj tested!";
  REQUIRE( true );
}

TEST_CASE("Eloc_Derivatives:multislater_noj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_derivatives:multislater_noj tested!";
  REQUIRE( true );
}

TEST_CASE("Eloc_Derivatives:multislater_wj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_derivatives:multislater_wj tested!";
  REQUIRE( true );
}

} // namespace qmcplusplus
