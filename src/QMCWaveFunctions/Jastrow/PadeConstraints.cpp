//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  ////////////////////////////////////////
  //PadeConstraints definitions
  ////////////////////////////////////////
  PadeConstraints::PadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
    {
      JComponent.set(MULTIPLE);
      JComponent.set(ONEBODY);
      JComponent.set(TWOBODY);
    }
   
  PadeConstraints::~PadeConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool PadeConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("B"));
    if(vit == inVars.end()) return false; //disaster, need to abort
    ID=(*vit).second.first; B=(*vit).second.second;
    return true;
  }

  OrbitalBase* PadeConstraints::createTwoBody() {
    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(targetPtcl);
    if(IgnoreSpin) {
      app_log() << "  PadeConstraints::Adding Spin-independent Pade Two-Body Jastrow " << endl;
      FuncType *func=new FuncType(-0.5,B);
      for(int i=0; i<4; i++) {
        J2->addFunc(func);
      }
      FuncList.push_back(func);
    } else {
      app_log() << "  PadeConstraints::Adding Spin-dependent Pade Two-Body Jastrow " << endl;
      FuncType *funcUU=new FuncType(-0.25,B);
      FuncType *funcUD=new FuncType(-0.5,B);
      J2->addFunc(funcUU);
      J2->addFunc(funcUD);
      J2->addFunc(funcUD);
      J2->addFunc(funcUU);

      FuncList.push_back(funcUU);
      FuncList.push_back(funcUD);
    }
    app_log() << "  PadeConstraints:: B = " << B <<endl;
    return J2;
  }

  OrbitalBase* PadeConstraints::createOneBody(ParticleSet& source) {
    app_log() << "  PadeBuilder::Adding Pade One-Body Jastrow with effective ionic charges." << endl;
    typedef OneBodyJastrow<FuncType> JneType;
    JneType* J1 = new JneType(source,targetPtcl);
    SpeciesSet& Species(source.getSpeciesSet());
    int ng=Species.getTotalNum();
    int icharge = Species.addAttribute("charge");
    for(int ig=0; ig<ng; ig++) {
      RealType zeff=Species(icharge,ig);
      app_log() << "    " << Species.speciesName[ig] <<  " Zeff = " << zeff << endl;
      FuncType *func=new FuncType(-zeff,B,std::pow(2*zeff,0.25));
      J1->addFunc(ig,func);
      FuncList.push_back(func);
    }
    return J1;
  }

  ////////////////////////////////////////
  //ScaledPadeConstraints definitions
  ////////////////////////////////////////
  ScaledPadeConstraints::ScaledPadeConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
    {
      JComponent.set(TWOBODY);
    }
  ScaledPadeConstraints::~ScaledPadeConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool ScaledPadeConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator bit(inVars.find("B"));
    map<string,pair<string,RealType> >::iterator cit(inVars.find("C"));
    if(bit == inVars.end() || cit == inVars.end()) return false; 
    BID=(*bit).second.first; B=(*bit).second.second; 
    CID=(*cit).second.first; C=(*cit).second.second; 
    return true;
  }

  OrbitalBase* ScaledPadeConstraints::createTwoBody() {
    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(targetPtcl);
    if(IgnoreSpin) {
      app_log() << "  ScaledPadeConstraints::Adding Spin-independent Pade Two-Body Jastrow " << endl;
      FuncType *func=new FuncType(-0.5,B,C);
      for(int i=0; i<4; i++) {
        J2->addFunc(func);
      }
      FuncList.push_back(func); 
    } else {
      app_log() << "  ScaledPadeConstraints::Adding Spin-dependent Pade Two-Body Jastrow " << endl;
      FuncType *funcUU=new FuncType(-0.25,B,C);
      FuncType *funcUD=new FuncType(-0.5,B,C);
      J2->addFunc(funcUU);
      J2->addFunc(funcUD);
      J2->addFunc(funcUD);
      J2->addFunc(funcUU);

      FuncList.push_back(funcUU);
      FuncList.push_back(funcUD);
    }
    app_log() << "  ScaledPadeConstraints:: B = " << B << " C = " << C <<endl;
    return J2;
  }

  OrbitalBase* 
    ScaledPadeConstraints::createOneBody(ParticleSet& source) {
    //return 0 for now
    return 0;
  }

  //////////////////////////////////////////
  ////PadeOnGridConstraints definitions
  //////////////////////////////////////////
  //PadeOnGridConstraints::~PadeOnGridConstraints() {
  //  delete_iter(FuncList.begin(), FuncList.end());
  //}

  //bool PadeOnGridConstraints::put(xmlNodePtr cur) {
  //  bool success=getVariables(cur);
  //  map<string,pair<string,RealType> >::iterator vit(inVars.find("B"));
  //  if(vit == inVars.end()) return false; //disaster, need to abort
  //  ID=(*vit).second.first; B=(*vit).second.second;
  //  return true;
  //}

  //void PadeOnGridConstraints::apply() {
  //  for(int i=0; i<FuncList.size(); i++) {
  //    InFuncList[i]->B0=B;
  //    FuncList[i]->reset();
  //  }
  //}

  //void PadeOnGridConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
  //  //potentially add Rcut
  //  outVars.add(ID,&B,1);
  //}

  //OrbitalBase* PadeOnGridConstraints::createTwoBody() {

  //  setRadialGrid(targetPtcl);

  //  typedef FuncType::FNOUT OutFuncType;
  //  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
  //  JeeType *J2 = new JeeType(target);
  //  if(IgnoreSpin) {
  //    app_log() << "  PadeOnGridConstraints::Adding Spin-independent Pade Two-Body Jastrow B=" << B << endl;
  //    //create an analytic input functor
  //    InFuncType *infunc=new InFuncType(-0.5,B);
  //    //create a numerical functor
  //    FuncType* nfunc= new FuncType;
  //    //initialize the numerical functor
  //    nfunc->initialize(infunc,myGrid,Rcut);

  //    InFuncList.push_back(infunc);
  //    FuncList.push_back(nfunc);
  //    for(int i=0; i<4; i++) J2->addFunc(nfunc);
  //  } else {
  //    app_log() << "  PadeOnGridConstraints::Adding Spin-dependent Pade Two-Body Jastrow B= " << B << endl;
  //    InFuncType *uu=new InFuncType(-0.25,B);
  //    InFuncType *ud=new InFuncType(-0.5,B);

  //    FuncType *funcUU= new FuncType; 
  //    funcUU->initialize(uu,myGrid,Rcut);

  //    FuncType *funcUD= new FuncType; 
  //    funcUD->initialize(ud,myGrid,Rcut);

  //    InFuncList.push_back(uu);
  //    InFuncList.push_back(ud);

  //    FuncList.push_back(funcUU);
  //    FuncList.push_back(funcUD);

  //    J2->addFunc(funcUU);//uu
  //    J2->addFunc(funcUD);//ud
  //    J2->addFunc(funcUD);//du
  //    J2->addFunc(funcUU);//dd
  //  }
  //  return J2;
  //}

  //OrbitalBase* PadeOnGridConstraints::createOneBody(ParticleSet& source) {
  //  //return 0 for the moment
  //  return 0;
  //}

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
