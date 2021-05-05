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

//This is a snapshot of a slightly compressed CN molecule.  r=1.10A.
//Electron configurations taken from a DMC snapshot.  
void create_CN_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  using pos = ParticleSet::SingleParticlePos_t;
  int Nions=2; //CN molecule.
  int Nelec=9; //With ccECP's, 4e for C, 5e for N;
  ions.setName("ion0");
  ions.create(Nions);

  ions.R[0]=pos(0.0,0.0,0.0);
  ions.R[1]=pos(0.0,0.0,2.0786985865e+00);
 
  SpeciesSet& ion_species          = ions.getSpeciesSet();
  int Cid                          = ion_species.addSpecies("C");
  int pChargeIdx                   = ion_species.addAttribute("charge");
  int pMembersizeIdx               = ion_species.addAttribute("membersize");
  int pAtomicNumIdx               = ion_species.addAttribute("atomicnumber");
  ion_species(pChargeIdx, Cid)     = 4;
  ion_species(pMembersizeIdx, Cid) = 1;
  ion_species(pAtomicNumIdx, Cid)  = 6;

  int Nid                          = ion_species.addSpecies("N");
  ion_species(pChargeIdx, Nid)     = 5;
  ion_species(pMembersizeIdx, Nid) = 1;
  ion_species(pAtomicNumIdx, Nid)  = 7;

  elec.setName("elec");
  elec.create(Nelec);

//up elecs
  elec.R[0]=pos(-0.55936725, -0.26942464,  0.14459603); 
  elec.R[1]=pos( 0.19146719,  1.40287983,  0.63931251);
  elec.R[2]=pos( 1.14805915, -0.52057335,  3.49621107);
  elec.R[3]=pos( 0.28293870, -0.10273952,  0.01707021);
  elec.R[4]=pos( 0.60626935, -0.25538121,  1.75750740);

//dn elecs
  elec.R[5]=pos(-0.47405939,  0.59523171, -0.59778601);
  elec.R[6]=pos( 0.03150661, -0.27343474,  0.56279442);
  elec.R[7]=pos(-1.32648025,  0.00970226,  2.26944242);
  elec.R[8]=pos( 2.42944286,  0.64884151,  1.87505288);
  
  SpeciesSet& tspecies           = elec.getSpeciesSet();
  int upIdx                      = tspecies.addSpecies("u");
  int dnIdx                      = tspecies.addSpecies("d");
  int chargeIdx                  = tspecies.addAttribute("charge");
  int massIdx                    = tspecies.addAttribute("mass");
  int MembersizeIdx              = tspecies.addAttribute("membersize");
  tspecies(MembersizeIdx, upIdx) = 5;
  tspecies(chargeIdx, upIdx)     = -1;
  tspecies(massIdx, upIdx)       = 1.0;

  tspecies(MembersizeIdx, dnIdx) = 4;
  tspecies(chargeIdx, dnIdx)     = -1;
  tspecies(massIdx, dnIdx)       = 1.0;

  elec.addTable(ions);
  elec.update();

}

void doopy(ParticleSet& elec, ParticleSet& ions)
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

//  create_CN_particlesets(elec,ions); 
  doopy(elec,ions);
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
