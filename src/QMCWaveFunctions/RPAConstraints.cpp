//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#include "QMCWaveFunctions/RPAConstraints.h"
#include "QMCWaveFunctions/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/OneBodyJastrowFunction.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  ////////////////////////////////////////
  //RPAConstraints definitions
  ////////////////////////////////////////
  RPAConstraints::~RPAConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool RPAConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID="rs"; Rs=-1.0;
        return true;//assign the default value
      }
    }
    ID=(*vit).second.first; Rs=(*vit).second.second;
    return true;
  }

  OrbitalBase* RPAConstraints::createTwoBody(ParticleSet& target) {

    if(Rs<0) {
      if(target.Lattice.BoxBConds[0]) {
        Rs=std::pow(3.0/4.0/M_PI*target.Lattice.Volume/static_cast<RealType>(target.getTotalNum()),1.0/3.0);
      } else {
        Rs=1.0;
      }
    }

    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(target);
    if(IgnoreSpin) {
      FuncType *func=new FuncType(false);
      func->reset(Rs);
      for(int i=0; i<4; i++) {
        J2->addFunc(func);
      }
      FuncList.push_back(func);
      app_log() << "  RPAConstraints::Adding Spin-independent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F=" << 1.0/func->B << endl;
    } else {
      FuncType *funcUU=new FuncType(true);
      FuncType *funcUD=new FuncType(false);
      funcUU->reset(Rs);
      funcUD->reset(Rs);

      J2->addFunc(funcUU);
      J2->addFunc(funcUD);
      J2->addFunc(funcUD);
      J2->addFunc(funcUU);

      FuncList.push_back(funcUU);
      FuncList.push_back(funcUD);
      app_log() << "  RPAConstraints::Adding Spin-dependent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F(uu)=" << 1.0/funcUU->B << "  F(ud)= " << 1.0/funcUD->B << endl;
    }
    return J2;
  }

  OrbitalBase* RPAConstraints::createOneBody(ParticleSet& target, ParticleSet& source) {
    return 0;
  }

  ////////////////////////////////////////
  //RPAPBCConstraints definitions
  ////////////////////////////////////////
  RPAPBCConstraints::~RPAPBCConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool RPAPBCConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID="rs"; Rs=-1.0;
        return true;//assign the default value
      }
    }
    ID=(*vit).second.first; Rs=(*vit).second.second;
    return true;
  }

  void RPAPBCConstraints::apply() {
    for(int i=0; i<FuncList.size(); i++) {
      InFuncList[i]->reset(Rs);
      FuncList[i]->reset();
    }
  }

  void RPAPBCConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
    //potentially add Rcut
    outVars.add(ID,&Rs,1);
  }

  OrbitalBase* RPAPBCConstraints::createTwoBody(ParticleSet& target) {

    setRadialGrid(target);

    if(Rs<0) {
      if(target.Lattice.BoxBConds[0]) {
        Rs=std::pow(3.0/4.0/M_PI*target.Lattice.Volume/static_cast<RealType>(target.getTotalNum()),1.0/3.0);
      } else {
        Rs=1.0;
      }
    }


    typedef FuncType::FNOUT OutFuncType;
    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(target);
    if(IgnoreSpin) {
      //create an analytic input functor
      InFuncType *infunc=new InFuncType(false);
      infunc->reset(Rs);
      //create a numerical functor
      FuncType* nfunc= new FuncType;
      //initialize the numerical functor
      nfunc->initialize(infunc,myGrid,Rcut);

      InFuncList.push_back(infunc);
      FuncList.push_back(nfunc);
      for(int i=0; i<4; i++) J2->addFunc(nfunc);

      app_log() << "  RPAPBCConstraints::Adding Spin-independent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F=" << 1.0/infunc->B << endl;
    } else {
      InFuncType *uu=new InFuncType(true);
      InFuncType *ud=new InFuncType(false);
      uu->reset(Rs);
      ud->reset(Rs);

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

      app_log() << "  RPAPBCConstraints::Adding Spin-dependent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F(uu)=" << 1.0/uu->B << "  F(ud)= " << 1.0/ud->B << endl;
    }
    return J2;
  }

  OrbitalBase* RPAPBCConstraints::createOneBody(ParticleSet& target, ParticleSet& source) {
    return 0;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
