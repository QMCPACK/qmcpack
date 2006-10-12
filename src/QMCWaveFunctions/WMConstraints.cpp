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
#include "QMCWaveFunctions/WMConstraints.h"
#include "QMCWaveFunctions/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/OneBodyJastrowFunction.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  WMConstraints::~WMConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool WMConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("B"));
    if(vit == inVars.end()) return false; //disaster, need to abort
    ID=(*vit).second.first; B=(*vit).second.second;
    return true;
  }

  void WMConstraints::apply() {
    for(int i=0; i<FuncList.size(); i++) {
      FuncList[i]->reset();
    }
  }

  void WMConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
    //outVars.add(ID,&B,1);
    //potentially add Rcut
  }

  OrbitalBase* WMConstraints::createTwoBody(ParticleSet& target) {

    setRadialGrid(target);

    typedef FuncType::FNOUT OutFuncType;
    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(target);
    if(IgnoreSpin) {
      app_log() << "  WMConstraints::Adding Spin-independent Pade Two-Body Jastrow B=" << B << endl;
      //create an analytic input functor
      InFuncType *infunc=new InFuncType(-0.5,B);
      //create a numerical functor
      FuncType* nfunc= new FuncType;
      //initialize the numerical functor
      nfunc->initialize(infunc,myGrid,Rcut);

      InFuncList.push_back(infunc);
      FuncList.push_back(nfunc);
      for(int i=0; i<4; i++) J2->addFunc(nfunc);
    } else {
      app_log() << "  WMConstraints::Adding Spin-dependent Pade Two-Body Jastrow B= " << B << endl;
      InFuncType *uu=new InFuncType(-0.25,B);
      InFuncType *ud=new InFuncType(-0.5,B);

      FuncType *funcUU= new FuncType; 
      funcUU->initialize(uu,myGrid,Rcut);

      FuncType *funcUD= new FuncType; 
      funcUD->initialize(ud,myGrid,Rcut);

      InFuncList.push_back(uu);
      InFuncList.push_back(ud);

      FuncList.push_back(funcUU);
      FuncList.push_back(funcUD);

      J2->addFunc(funcUU);//uu
      J2->addFunc(funcUD);//ud
      J2->addFunc(funcUD);//du
      J2->addFunc(funcUU);//dd
    }
    return J2;
  }

  OrbitalBase* WMConstraints::createOneBody(ParticleSet& target, ParticleSet& source) {
    //return 0 for the moment
    return 0;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
