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
#include "Utilities/RandomGenerator.h"
#include "QMCWaveFunctions/TWFPrototype.h"
#include "Numerics/DeterminantOperators.h"

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
  app_log() << "====Ion Derivative Test: Single Slater No Jastrow====\n";
  using RealType = QMCTraits::RealType;

  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec, ions);

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
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
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
  //HamTest   Total -1.617071104264908143e+01
  //HamTest Kinetic 9.182193792758450712e+00
  //HamTest ElecElec 1.901556057075800865e+01
  //HamTest IonIon 9.621404531608845900e+00
  //HamTest LocalECP -6.783942829945100073e+01
  //HamTest NonLocalECP 1.384955836167661225e+01


  RealType logpsi = psi->evaluateLog(elec);
  REQUIRE(logpsi == Approx(-14.233853149));

  QMCHamiltonian* ham = hf.getH();

  RealType eloc = ham->evaluateDeterministic(elec);
  enum observ_id
  {
    KINETIC = 0,
    ELECELEC,
    IONION,
    LOCALECP,
    NONLOCALECP
  };
  REQUIRE(eloc == Approx(-1.6170527168e+01));
  REQUIRE(ham->getObservable(ELECELEC) == Approx(1.9015560571e+01));
  REQUIRE(ham->getObservable(IONION) == Approx(9.6214045316e+00));
  REQUIRE(ham->getObservable(LOCALECP) == Approx(-6.7839428299e+01));
  REQUIRE(ham->getObservable(KINETIC) == Approx(9.1821937928e+00));
  REQUIRE(ham->getObservable(NONLOCALECP) == Approx(1.3849558361e+01));

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
  REQUIRE(wf_grad[0][0] == Approx(-1.9044650674260308));
  REQUIRE(wf_grad[0][1] == Approx(2.1257764985627148));
  REQUIRE(wf_grad[0][2] == Approx(7.0556319963444016));
  REQUIRE(wf_grad[1][0] == Approx(1.4233346821157509));
  REQUIRE(wf_grad[1][1] == Approx(-0.1446706081154048));
  REQUIRE(wf_grad[1][2] == Approx(0.1440176276013005));

  //Kinetic Force
  hf_term    = 0.0;
  pulay_term = 0.0;
  (ham->getHamiltonian(KINETIC))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
#if defined(MIXED_PRECISION)
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(1.0852823603357820).epsilon(1e-4));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(24.2154119471038562).epsilon(1e-4));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(111.8849364775797852).epsilon(1e-4));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.1572063443997536).epsilon(1e-4));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(-3.3743242489947529).epsilon(1e-4));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(7.5625192454964454).epsilon(1e-4));
#else
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(1.0852823603357820));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(24.2154119471038562));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(111.8849364775797852));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.1572063443997536));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(-3.3743242489947529));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(7.5625192454964454));
#endif
  //NLPP Force
  hf_term    = 0.0;
  pulay_term = 0.0;
  double val =
      (ham->getHamiltonian(NONLOCALECP))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);

//MP fails the first REQUIRE with 24.22544.  Just bypass the checks in those builds.
#if defined(MIXED_PRECISION)
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(24.2239540340527491).epsilon(2e-4));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(-41.9981344310649263).epsilon(2e-4));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(-98.9123955744908159).epsilon(2e-4));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.5105943834091704).epsilon(2e-4));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(1.1345766918857692).epsilon(2e-4));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(-5.2293234395150989).epsilon(2e-4));
#else
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(24.2239540340527491));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(-41.9981344310649263));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(-98.9123955744908159));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.5105943834091704));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(1.1345766918857692));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(-5.2293234395150989));
#endif

 //End of deterministic tests.  Let's call evaluateIonDerivs and evaluateIonDerivsDeterministic at the
 //QMCHamiltonian level to make sure there are no crashes.  
 
  hf_term    = 0.0;
  pulay_term = 0.0;
  wf_grad    = 0.0;
  ham->evaluateIonDerivsDeterministic(elec,ions,*psi,hf_term,pulay_term,wf_grad);
 
  REQUIRE(dot(hf_term[0],hf_term[0]) != Approx(0));
  REQUIRE(dot(pulay_term[0],pulay_term[0]) != Approx(0));
  REQUIRE(dot(wf_grad[0],wf_grad[0]) != Approx(0));
 
  REQUIRE(dot(hf_term[1],hf_term[1]) != Approx(0));
  REQUIRE(dot(pulay_term[1],pulay_term[1]) != Approx(0));
  REQUIRE(dot(wf_grad[1],wf_grad[1]) != Approx(0));
 
  hf_term    = 0.0;
  pulay_term = 0.0;
  wf_grad    = 0.0;
  RandomGenerator_t myrng;
  ham->setRandomGenerator(&myrng);
  ham->evaluateIonDerivs(elec,ions,*psi,hf_term,pulay_term,wf_grad);
  
  REQUIRE(dot(hf_term[0],hf_term[0]) != Approx(0));
  REQUIRE(dot(pulay_term[0],pulay_term[0]) != Approx(0));
  REQUIRE(dot(wf_grad[0],wf_grad[0]) != Approx(0));
 
  REQUIRE(dot(hf_term[1],hf_term[1]) != Approx(0));
  REQUIRE(dot(pulay_term[1],pulay_term[1]) != Approx(0));
  REQUIRE(dot(wf_grad[1],wf_grad[1]) != Approx(0));
 
}

TEST_CASE("Eloc_Derivatives:proto_sd_wj","[hamiltonian]")
{
  app_log() <<"========================================================================================\n";
  app_log() <<"========================================================================================\n";
  app_log() <<"====================Ion Derivative Test:  Prototype Single Slater+Jastrow===============\n";
  app_log() <<"========================================================================================\n";
  app_log() <<"========================================================================================\n";
  using RealType = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;  
  Communicate* c;
  c= OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec,ions);
    
  Libxml2Document doc2;
  bool okay = doc2.parse("cn.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  HamiltonianFactory::PsiPoolType psi_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  WaveFunctionFactory wff("psi0", elec, particle_set_map, c);
  psi_map["psi0"] = &wff;

  SPOSetBuilderFactory bf(c, elec, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);


  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//sposet", doc2.getXPathContext());
  SPOSet* sposet_up = bb->createSPOSet(slater_base[0]);
  SPOSet* sposet_dn = bb->createSPOSet(slater_base[1]);

    
  TWFPrototype twf;
  twf.add_determinant(elec,0,sposet_up);
  twf.add_determinant(elec,1,sposet_dn);

  SPOSet::ValueVector_t values;
  SPOSet::GradVector_t dpsi;
  SPOSet::ValueVector_t d2psi;
  values.resize(9);
  dpsi.resize(9);
  d2psi.resize(9);


  HamiltonianFactory hf("h0", elec, particle_set_map, psi_map, c);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
           <pseudo elementType=\"N\" href=\"N.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";

  Libxml2Document doc;
  okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  hf.put(root);

  QMCHamiltonian* ham = hf.getH();

  enum observ_id
  {
    KINETIC = 0,
    ELECELEC,
    IONION,
    LOCALECP,
    NONLOCALECP
  };

  using ValueMatrix_t = SPOSet::ValueMatrix_t;

  int IONINDEX=1;

  ValueMatrix_t upmat;
  ValueMatrix_t dnmat;
  upmat.resize(5,14);
  dnmat.resize(4,14);

  std::vector<ValueMatrix_t> matlist;
  std::vector<ValueMatrix_t> Bkin;
  std::vector<std::vector<ValueMatrix_t> > dM;
  std::vector<std::vector<ValueMatrix_t> > dB;
  matlist.push_back(upmat);
  matlist.push_back(dnmat);
  
  dM.push_back(matlist);
  dM.push_back(matlist);
  dM.push_back(matlist);

  dB.push_back(matlist);
  dB.push_back(matlist);
  dB.push_back(matlist);

  Bkin.push_back(upmat);
  Bkin.push_back(dnmat);

  twf.get_M(elec,matlist);
  
  OperatorBase* kinop=ham->getHamiltonian(KINETIC);
 
  kinop-> evaluateOneBodyOpMatrix(elec,twf,Bkin);

  std::vector<ValueMatrix_t> minv;
  std::vector<ValueMatrix_t> Bkin_gs;
  std::vector<std::vector<ValueMatrix_t> > dB_gs;
  std::vector<std::vector<ValueMatrix_t> > dM_gs;
  std::vector<ValueMatrix_t> tmp_gs;
  //This handles the computation of the ground state determinant.  Here for now.
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueMatrix_t gs_m;
    ValueMatrix_t gs_kin;
    ValueMatrix_t gs_minv;

    gs_m.resize(ptclnum,ptclnum);
    gs_kin.resize(ptclnum,ptclnum);
    gs_minv.resize(ptclnum,ptclnum);
    tmp_gs.push_back(gs_m); 
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
      {
        gs_m[i][j]=matlist[id][i][j];
        gs_kin[i][j]=Bkin[id][i][j];
      }
    gs_minv=gs_m;
    invert_matrix(gs_minv);
    minv.push_back(gs_minv);
    Bkin_gs.push_back(gs_kin);
  }
  dB_gs.push_back(tmp_gs);
  dB_gs.push_back(tmp_gs);
  dB_gs.push_back(tmp_gs);

  dM_gs.push_back(tmp_gs);
  dM_gs.push_back(tmp_gs);
  dM_gs.push_back(tmp_gs);

  kinop->evaluateOneBodyOpMatrixForceDeriv(elec,ions,twf,IONINDEX,dB);
  twf.get_igrad_M(elec,ions,IONINDEX,dM);
  for(int idim=0; idim<OHMMS_DIM; idim++)
    for(int id=0; id<matlist.size(); id++)
    {
      int ptclnum = twf.num_particles(id);
      for(int i=0; i<ptclnum; i++)
        for(int j=0; j<ptclnum; j++)
        {      
          dB_gs[idim][id][i][j] = dB[idim][id][i][j];
          dM_gs[idim][id][i][j] = dM[idim][id][i][j];
        }
    }
  ValueType keval = 0.0;
  //Now to compute the kinetic energy
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueType ke_id = 0.0;
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
      {
        ke_id+=minv[id][i][j]*Bkin[id][j][i];
      }
    keval+=ke_id;
  }

  app_log()<<" KEVal = "<<keval<<std::endl;
 
  //Now we do the forces.  
  //First build up the intermediate matrices.
  std::vector<ValueMatrix_t> X;
    
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueMatrix_t Xid;
    ValueMatrix_t tmpmat;

    Xid.resize(ptclnum,ptclnum);
    tmpmat.resize(ptclnum,ptclnum);
    //(B*A^-1)
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
        for(int k=0; k<ptclnum; k++)
        {
          tmpmat[i][j]+=Bkin_gs[id][i][k]*minv[id][k][j];
        } 

    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
        for(int k=0; k<ptclnum; k++)
        {
          Xid[i][j]+=minv[id][i][k]*tmpmat[k][j];
        } 
    X.push_back(Xid);
  }

  for(int idim=0; idim<OHMMS_DIM; idim++)
  {
    ValueType dval=0.0;
    ValueType dwfn=0.0;
    for(int id=0; id<matlist.size(); id++)
    {
      int ptclnum = twf.num_particles(id);
      ValueType dval_id=0.0; 
      ValueType dwfn_id=0.0;
      for(int i=0; i<ptclnum; i++)
        for(int j=0; j<ptclnum; j++)
        {
          dval_id+=minv[id][i][j]*dB_gs[idim][id][j][i]-X[id][i][j]*dM_gs[idim][id][j][i];
          dwfn_id+=minv[id][i][j]*dM_gs[idim][id][j][i];
        }
      dval+=dval_id;
      dwfn+=dwfn_id;
    }
    app_log()<<"F["<<IONINDEX<<"]["<<idim<<"]="<<dval<<std::endl;
    app_log()<<"dTWF["<<IONINDEX<<"]["<<idim<<"]="<<dwfn<<std::endl;
  } 

  app_log()<<" Now evaluating nonlocalecp\n"; 
  OperatorBase* nlppop=ham->getHamiltonian(NONLOCALECP);
  app_log()<<"nlppop = "<<nlppop<<std::endl;
  app_log()<<"  Evaluated.  Calling evaluteOneBodyOpMatrix\n";
  
  for(int id; id<matlist.size(); id++)
    Bkin[id]=0.0;

  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bkin);
  
  app_log()<<" Done\n";
  //This handles the computation of the ground state determinant.  Here for now.
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueMatrix_t gs_m;
    ValueMatrix_t gs_kin;
    ValueMatrix_t gs_minv;

    gs_m.resize(ptclnum,ptclnum);
    gs_kin.resize(ptclnum,ptclnum);
    gs_minv.resize(ptclnum,ptclnum);
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
      {
        gs_m[i][j]=matlist[id][i][j];
        gs_kin[i][j]=Bkin[id][i][j];
      }
    gs_minv=gs_m;
    invert_matrix(gs_minv);
    minv.push_back(gs_minv);
    Bkin_gs.push_back(gs_kin);
  }
  
  ValueType nlpp = 0.0;
  //Now to compute the kinetic energy
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueType nlpp_id = 0.0;
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
      {
        nlpp_id+=minv[id][i][j]*Bkin[id][j][i];
      }
    nlpp+=nlpp_id;
  }
   
  app_log()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log()<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log()<<" NLPPVal = "<<nlpp<<std::endl;

  std::vector<ValueMatrix_t> Bplus;
  std::vector<ValueMatrix_t> Bminus;
  
  for(int i=0; i<matlist.size(); i++)
  {
    Bplus.push_back(Bkin[i]);
    Bminus.push_back(Bkin[i]);
 //   app_log()<<" Species "<<i<<std::endl;
//    app_log()<<Bkin[i]<<std::endl<<std::endl;
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  for(int i=0; i<matlist.size(); i++)
  {
    Bplus[i]=0.0;
    Bminus[i]=0.0;
  }

  RealType delta=1e-4;
 

  RealType xold=ions.R[IONINDEX][0];
  ions.R[IONINDEX][0]=xold+delta;

  ions.update();
  elec.update();

  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bplus);
  ions.R[IONINDEX][0]=xold-delta;

  ions.update();
  elec.update();
  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bminus);
  {std::vector<ValueMatrix_t> Bdiff_vec; 
  for(int i=0; i<matlist.size(); i++)
  {
    ValueMatrix_t Bdiff;
    int ptclnum = twf.num_particles(i);
    int norbs = twf.num_orbitals(i);
    Bdiff.resize(ptclnum,norbs);
    Bdiff=(Bplus[i]-Bminus[i])*(1.0/2.0/delta);
    Bdiff_vec.push_back(Bdiff);
    app_log()<<" [0] Species BDIFF "<<i<<std::endl;
    app_log()<<Bdiff<<std::endl;
  }}
  
  ions.R[IONINDEX][0]=xold;
  ions.update();
  elec.update();

  for(int i=0; i<matlist.size(); i++)
  {
    Bplus[i]=0.0;
    Bminus[i]=0.0;
  }

 
  xold=ions.R[IONINDEX][1];
  ions.R[IONINDEX][1]=xold+delta;

  ions.update();
  elec.update();
  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bplus);
  ions.R[IONINDEX][1]=xold-delta;
  ions.update();
  elec.update();

  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bminus);
  {std::vector<ValueMatrix_t> Bdiff_vec; 
  for(int i=0; i<matlist.size(); i++)
  {
    ValueMatrix_t Bdiff;
    int ptclnum = twf.num_particles(i);
    int norbs = twf.num_orbitals(i);
    Bdiff.resize(ptclnum,norbs);
    Bdiff=(Bplus[i]-Bminus[i])*(1.0/2.0/delta);
    Bdiff_vec.push_back(Bdiff);
    app_log()<<"[1] Species BDIFF "<<i<<std::endl;
    app_log()<<Bdiff<<std::endl;
  }}
  
  ions.R[IONINDEX][1]=xold;
  ions.update();
  elec.update();
  
  xold=ions.R[IONINDEX][2];
  ions.R[IONINDEX][2]=xold+delta;

  ions.update();
  elec.update();

  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bplus);
  ions.R[0][2]=xold-delta;
  ions.update();
  elec.update();

  nlppop-> evaluateOneBodyOpMatrix(elec,twf,Bminus);
  {std::vector<ValueMatrix_t> Bdiff_vec; 
  for(int i=0; i<matlist.size(); i++)
  {
    ValueMatrix_t Bdiff;
    int ptclnum = twf.num_particles(i);
    int norbs = twf.num_orbitals(i);
    Bdiff.resize(ptclnum,norbs);
    Bdiff=(Bplus[i]-Bminus[i])*(1.0/2.0/delta);
    Bdiff_vec.push_back(Bdiff);
    app_log()<<"[2] Species BDIFF "<<i<<std::endl;
    app_log()<<Bdiff<<std::endl;
  }}
   
  ions.R[IONINDEX][2]=xold;
  ions.update();
  elec.update();

 */ 
  for(int idim=0; idim<OHMMS_DIM; idim++)
    for(int id=0; id<matlist.size(); id++)
      dB[idim][id]=0.0;
  
  nlppop->evaluateOneBodyOpMatrixForceDeriv(elec,ions,twf,IONINDEX,dB);
  for(int idim=0; idim<OHMMS_DIM; idim++)
    for(int id=0; id<matlist.size(); id++)
    {
      int ptclnum = twf.num_particles(id);
      for(int i=0; i<ptclnum; i++)
        for(int j=0; j<ptclnum; j++)
        {      
          dB_gs[idim][id][i][j] = dB[idim][id][i][j];
          //dB_gs[0][id][i][j] = Bdiff_vec[id][i][j];
          Bkin_gs[id][i][j] = Bkin[id][i][j];
        }
    }
  
  std::vector<ValueMatrix_t> X2;
  for(int id=0; id<matlist.size(); id++)
  {
    int ptclnum = twf.num_particles(id);
    ValueMatrix_t Xid;
    ValueMatrix_t tmpmat;
    Xid.resize(ptclnum,ptclnum);
    Xid=0.0;
    tmpmat=0.0;
    tmpmat.resize(ptclnum,ptclnum);
    //(B*A^-1)
    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
        for(int k=0; k<ptclnum; k++)
        {
          tmpmat[i][j]+=Bkin_gs[id][i][k]*minv[id][k][j];
        } 

    for(int i=0; i<ptclnum; i++)
      for(int j=0; j<ptclnum; j++)
        for(int k=0; k<ptclnum; k++)
        {
          Xid[i][j]+=minv[id][i][k]*tmpmat[k][j];
        } 
    X2.push_back(Xid);
  }

 // int idim=0;
  for(int idim=0; idim<OHMMS_DIM; idim++)
  {
    ValueType dval=0.0;
    ValueType dwfn=0.0;
    for(int id=0; id<matlist.size(); id++)
    {
    //  app_log()<<" BDiff Species "<<id<<std::endl;
    //  app_log()<<dB[idim][id]<<std::endl<<std::endl;
      int ptclnum = twf.num_particles(id);
      ValueType dval_id=0.0; 
      ValueType dwfn_id=0.0;
      for(int i=0; i<ptclnum; i++)
        for(int j=0; j<ptclnum; j++)
        {
          dval_id+=minv[id][i][j]*dB_gs[idim][id][j][i]-X2[id][i][j]*dM_gs[idim][id][j][i];
          dwfn_id+=minv[id][i][j]*dM_gs[idim][id][j][i];
        }
      dval+=dval_id;
      dwfn+=dwfn_id;
    }
    app_log()<<"F["<<IONINDEX<<"]["<<idim<<"]="<<dval<<std::endl;
    app_log()<<"dTWF["<<IONINDEX<<"]["<<idim<<"]="<<dwfn<<std::endl;
  } 

   
}

/*TEST_CASE("Eloc_Derivatives:slater_wj", "[hamiltonian]")
{
  app_log() << "====Ion Derivative Test: Single Slater+Jastrow====\n";
  using RealType = QMCTraits::RealType;

  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec, ions);

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
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
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
  bool wfokay = wfdoc.parse("cn.wfj.xml");
  REQUIRE(wfokay);

  xmlNodePtr wfroot = wfdoc.getRoot();
  wff.put(wfroot);

  TrialWaveFunction* psi = wff.getTWF();
  REQUIRE(psi != nullptr);

  //Output of WFTester Eloc test for this ion/electron configuration.
  //  Logpsi: (-8.945509461103977600e+00,0.000000000000000000e+00)
  //  HamTest   Total -1.779268125690864721e+01
  //  HamTest Kinetic 7.673240415372170276e+00
  //  HamTest ElecElec 1.901556057075800865e+01
  //  HamTest IonIon 9.621404531608845900e+00
  //  HamTest LocalECP -6.783942829945100073e+01
  //  HamTest NonLocalECP 1.373654152480333224e+01

  RealType logpsi = psi->evaluateLog(elec);
  REQUIRE(logpsi == Approx(-8.9455094611e+00));

  QMCHamiltonian* ham = hf.getH();

  RealType eloc = ham->evaluateDeterministic(elec);
  enum observ_id
  {
    KINETIC = 0,
    ELECELEC,
    IONION,
    LOCALECP,
    NONLOCALECP
  };
  REQUIRE(eloc == Approx(-1.77926812569e+01));
  REQUIRE(ham->getObservable(ELECELEC) == Approx(1.9015560571e+01));
  REQUIRE(ham->getObservable(IONION) == Approx(9.6214045316e+00));
  REQUIRE(ham->getObservable(LOCALECP) == Approx(-6.7839428299e+01));
  REQUIRE(ham->getObservable(KINETIC) == Approx(7.6732404154e+00));
  REQUIRE(ham->getObservable(NONLOCALECP) == Approx(1.37365415248e+01));

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
  REQUIRE(wf_grad[0][0] == Approx(-1.8996878390353797));
  REQUIRE(wf_grad[0][1] == Approx(2.3247646590007776));
  REQUIRE(wf_grad[0][2] == Approx(7.9587196049502031));
  REQUIRE(wf_grad[1][0] == Approx(1.8093817104158914));
  REQUIRE(wf_grad[1][1] == Approx(-0.0966225639942308));
  REQUIRE(wf_grad[1][2] == Approx(-1.5197874544625731));

  //Kinetic Force
  hf_term    = 0.0;
  pulay_term = 0.0;
  (ham->getHamiltonian(KINETIC))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
#if defined(MIXED_PRECISION)
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(-3.3359153349010735).epsilon(1e-4));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(30.0487085581835309).epsilon(1e-4));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(126.5885230360197369).epsilon(1e-4));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.7271604366774223).epsilon(1e-4));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(-3.5321234918228579).epsilon(1e-4));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(5.8844148870917925).epsilon(1e-4));
#else
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(-3.3359153349010735));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(30.0487085581835309));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(126.5885230360197369));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.7271604366774223));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(-3.5321234918228579));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(5.8844148870917925));
#endif
  //NLPP Force
  hf_term    = 0.0;
  pulay_term = 0.0;
  double val =
      (ham->getHamiltonian(NONLOCALECP))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
//MP fails the first REQUIRE with 27.15313.  Just bypass the checks in those builds.
#if defined(MIXED_PRECISION)
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(27.1517161490208956).epsilon(2e-4));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(-42.8268964286715459).epsilon(2e-4));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(-101.5046844660360961).epsilon(2e-4));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.2255825024686260).epsilon(2e-4));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(1.1362118534918864).epsilon(2e-4));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(-4.5825638607333019).epsilon(2e-4));
#else
  REQUIRE(hf_term[0][0] + pulay_term[0][0] == Approx(27.1517161490208956));
  REQUIRE(hf_term[0][1] + pulay_term[0][1] == Approx(-42.8268964286715459));
  REQUIRE(hf_term[0][2] + pulay_term[0][2] == Approx(-101.5046844660360961));
  REQUIRE(hf_term[1][0] + pulay_term[1][0] == Approx(2.2255825024686260));
  REQUIRE(hf_term[1][1] + pulay_term[1][1] == Approx(1.1362118534918864));
  REQUIRE(hf_term[1][2] + pulay_term[1][2] == Approx(-4.5825638607333019));
#endif

 //End of deterministic tests.  Let's call evaluateIonDerivs and evaluateIonDerivsDeterministic at the
 //QMCHamiltonian level to make sure there are no crashes.  
 
  hf_term    = 0.0;
  pulay_term = 0.0;
  wf_grad    = 0.0;
  ham->evaluateIonDerivsDeterministic(elec,ions,*psi,hf_term,pulay_term,wf_grad);
 
  REQUIRE(dot(hf_term[0],hf_term[0]) != Approx(0));
  REQUIRE(dot(pulay_term[0],pulay_term[0]) != Approx(0));
  REQUIRE(dot(wf_grad[0],wf_grad[0]) != Approx(0));
 
  REQUIRE(dot(hf_term[1],hf_term[1]) != Approx(0));
  REQUIRE(dot(pulay_term[1],pulay_term[1]) != Approx(0));
  REQUIRE(dot(wf_grad[1],wf_grad[1]) != Approx(0));
 
  hf_term    = 0.0;
  pulay_term = 0.0;
  wf_grad    = 0.0;
  RandomGenerator_t myrng;
  ham->setRandomGenerator(&myrng);
  ham->evaluateIonDerivs(elec,ions,*psi,hf_term,pulay_term,wf_grad);
  
  REQUIRE(dot(hf_term[0],hf_term[0]) != Approx(0));
  REQUIRE(dot(pulay_term[0],pulay_term[0]) != Approx(0));
  REQUIRE(dot(wf_grad[0],wf_grad[0]) != Approx(0));
 
  REQUIRE(dot(hf_term[1],hf_term[1]) != Approx(0));
  REQUIRE(dot(pulay_term[1],pulay_term[1]) != Approx(0));
  REQUIRE(dot(wf_grad[1],wf_grad[1]) != Approx(0));
 
}*/

/*TEST_CASE("Eloc_Derivatives:multislater_noj", "[hamiltonian]")
{
  app_log() << "====Ion Derivative Test: Multislater No Jastrow====\n";
  using RealType = QMCTraits::RealType;

  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec, ions);

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
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
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
  bool wfokay = wfdoc.parse("cn.msd-wfnoj.xml");
  REQUIRE(wfokay);

  xmlNodePtr wfroot = wfdoc.getRoot();
  wff.put(wfroot);

  TrialWaveFunction* psi = wff.getTWF();
  REQUIRE(psi != nullptr);

  //Output of WFTester Eloc test for this ion/electron configuration.
  //Logpsi: (-1.411499619826623686e+01,0.000000000000000000e+00)
  //HamTest   Total -1.597690575658561407e+01
  //HamTest Kinetic 1.053500867576629574e+01
  //HamTest ElecElec 1.901556057075800865e+01
  //HamTest IonIon 9.621404531608845900e+00
  //HamTest LocalECP -6.783942829945100073e+01
  //HamTest NonLocalECP 1.269054876473223636e+01

  RealType logpsi = psi->evaluateLog(elec);
  REQUIRE(logpsi == Approx(-1.41149961982e+01));

  QMCHamiltonian* ham = hf.getH();

  RealType eloc = ham->evaluateDeterministic(elec);
  enum observ_id
  {
    KINETIC = 0,
    ELECELEC,
    IONION,
    LOCALECP,
    NONLOCALECP
  };
  REQUIRE(eloc == Approx(-1.59769057565e+01));
  REQUIRE(ham->getObservable(ELECELEC) == Approx(1.90155605707e+01));
  REQUIRE(ham->getObservable(IONION) == Approx(9.62140453161e+00));
  REQUIRE(ham->getObservable(LOCALECP) == Approx(-6.78394282995e+01));
  REQUIRE(ham->getObservable(KINETIC) == Approx(1.05350086757e+01));
  REQUIRE(ham->getObservable(NONLOCALECP) == Approx(1.26905487647e+01));

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

  //This is not implemented yet.  Uncomment to perform check after implementation.
  //Reference from finite differences on this configuration.
  //REQUIRE( wf_grad[0][0] == Approx(-1.7045200053189544));
  //REQUIRE( wf_grad[0][1] == Approx( 2.6980932676501368)); 
  //REQUIRE( wf_grad[0][2] == Approx( 6.5358393587011667)); 
  //REQUIRE( wf_grad[1][0] == Approx( 1.6322817486980055)); 
  //REQUIRE( wf_grad[1][1] == Approx( 0.0091648450606385));
  //REQUIRE( wf_grad[1][2] == Approx( 0.1031883398283639)); 


  //This is not implemented yet.  Uncomment to perform check after implementation.
  //Kinetic Force
  //hf_term=0.0;
  //pulay_term=0.0;
  //(ham->getHamiltonian(KINETIC))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
//#if defined(MIXED_PRECISION) 
  //REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 7.4631825180304636).epsilon(1e-4));
  //REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx( 26.0975954772035799).epsilon(1e-4));
  //REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx( 90.1646424427582218).epsilon(1e-4));
  //REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx( 3.8414153131327562).epsilon(1e-4));
  //REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(-2.3504392874684754).epsilon(1e-4));
  //REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx( 4.7454048248241065) .epsilon(1e-4));
//#else
  //REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 7.4631825180304636));
  //REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx( 26.0975954772035799));
  //REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx( 90.1646424427582218));
  //REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx( 3.8414153131327562));
  //REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(-2.3504392874684754));
  //REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx( 4.7454048248241065)); 
//#endif 
  //This is not implemented yet.  Uncomment to perform check after implementation.
  //NLPP Force
 // hf_term=0.0;
//  pulay_term=0.0;
//  double val=(ham->getHamiltonian(NONLOCALECP))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
//#if defined(MIXED_PRECISION)
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 18.9414437404167302).epsilon(2e-4));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(-42.9017371899931277).epsilon(2e-4));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(-78.3304792483008328).epsilon(2e-4));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx( 1.2122162598160457).epsilon(2e-4));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(-0.6163169101291999).epsilon(2e-4));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(-3.2996553033015625).epsilon(2e-4));
//#else
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 18.9414437404167302));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(-42.9017371899931277));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(-78.3304792483008328));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx( 1.2122162598160457));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx(-0.6163169101291999));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(-3.2996553033015625));
//#endif 
}*/

/*TEST_CASE("Eloc_Derivatives:multislater_wj", "[hamiltonian]")
{
  app_log() << "====Ion Derivative Test: Multislater+Jastrow====\n";
  using RealType = QMCTraits::RealType;

  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions;
  ParticleSet elec;

  create_CN_particlesets(elec, ions);

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
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
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
  bool wfokay = wfdoc.parse("cn.msd-wfj.xml");
  REQUIRE(wfokay);

  xmlNodePtr wfroot = wfdoc.getRoot();
  wff.put(wfroot);

  TrialWaveFunction* psi = wff.getTWF();
  REQUIRE(psi != nullptr);

  //Output of WFTester Eloc test for this ion/electron configuration.
  //Logpsi: (-8.693299948465634586e+00,0.000000000000000000e+00)
  //HamTest   Total -1.752111246795173116e+01
  //HamTest Kinetic 9.187721666379577101e+00
  //HamTest ElecElec 1.901556057075800865e+01
  //HamTest IonIon 9.621404531608845900e+00
  //HamTest LocalECP -6.783942829945100073e+01
  //HamTest NonLocalECP 1.249362906275283969e+01


  RealType logpsi = psi->evaluateLog(elec);
  REQUIRE(logpsi == Approx(-8.69329994846e+00));

  QMCHamiltonian* ham = hf.getH();

  RealType eloc = ham->evaluateDeterministic(elec);
  enum observ_id
  {
    KINETIC = 0,
    ELECELEC,
    IONION,
    LOCALECP,
    NONLOCALECP
  };
  REQUIRE(eloc == Approx(-1.75211124679e+01));
  REQUIRE(ham->getObservable(ELECELEC) == Approx(1.90155605707e+01));
  REQUIRE(ham->getObservable(IONION) == Approx(9.62140453161e+00));
  REQUIRE(ham->getObservable(LOCALECP) == Approx(-6.78394282995e+01));
  REQUIRE(ham->getObservable(KINETIC) == Approx(9.18772166638e+00));
  REQUIRE(ham->getObservable(NONLOCALECP) == Approx(1.24936290628e+01));

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

  //This is not implemented yet.  Uncomment to perform check after implementation.
  //Reference from finite differences on this configuration.
//  REQUIRE( wf_grad[0][0] == Approx(-1.7052805961093040));
//  REQUIRE( wf_grad[0][1] == Approx( 2.8914116872336133)); 
//  REQUIRE( wf_grad[0][2] == Approx( 7.3963610874194776)); 
//  REQUIRE( wf_grad[1][0] == Approx( 2.0450537814298286)); 
//  REQUIRE( wf_grad[1][1] == Approx( 0.0742023428479399));
//  REQUIRE( wf_grad[1][2] == Approx(-1.6411356565271260));

  //This is not implemented yet.  Uncomment to perform check after implementation.
  //Kinetic Force
//  hf_term=0.0;
//  pulay_term=0.0;
//  (ham->getHamiltonian(KINETIC))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
//#if defined(MIXED_PRECISION)
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx(  4.1783687883878429).epsilon(1e-4));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx( 32.2193450745800192).epsilon(1e-4));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(102.0214857307521896).epsilon(1e-4));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(  4.5063296809644271).epsilon(1e-4));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx( -2.3360060461996568).epsilon(1e-4));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(  2.9502526588842666).epsilon(1e-4));
//#else
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx(  4.1783687883878429));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx( 32.2193450745800192));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(102.0214857307521896));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(  4.5063296809644271));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx( -2.3360060461996568));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx(  2.9502526588842666)); 
//#endif 
  //This is not implemented yet.  Uncomment to perform check after implementation.
  //NLPP Force
//  hf_term=0.0;
//  pulay_term=0.0;
//  double val=(ham->getHamiltonian(NONLOCALECP))->evaluateWithIonDerivsDeterministic(elec, ions, *psi, hf_term, pulay_term);
//#if defined(MIXED_PRECISION)
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 21.6829856774403140).epsilon(2e-4));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(-43.4432406419382673).epsilon(2e-4));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(-80.1356331911584618).epsilon(2e-4));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(  0.9915030925178313).epsilon(2e-4));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx( -0.6012127592214256).epsilon(2e-4));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx( -2.7937129314814508).epsilon(2e-4));
//#else
//  REQUIRE( hf_term[0][0]+pulay_term[0][0] == Approx( 21.6829856774403140));
//  REQUIRE( hf_term[0][1]+pulay_term[0][1] == Approx(-43.4432406419382673));
//  REQUIRE( hf_term[0][2]+pulay_term[0][2] == Approx(-80.1356331911584618));
//  REQUIRE( hf_term[1][0]+pulay_term[1][0] == Approx(  0.9915030925178313));
//  REQUIRE( hf_term[1][1]+pulay_term[1][1] == Approx( -0.6012127592214256));
//  REQUIRE( hf_term[1][2]+pulay_term[1][2] == Approx( -2.7937129314814508)); 
//#endif 
}*/

} // namespace qmcplusplus
