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
#include "QMCWaveFunctions/LRTwoBodyJastrow.h"
#include "Utilities/IteratorUtility.h"
#include "LongRange/LRJastrowSingleton.h"

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
  RPAPBCConstraints::~RPAPBCConstraints() { ; }

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

  void RPAPBCConstraints::apply() { ; }

  void RPAPBCConstraints::addOptimizables(VarRegistry<RealType>& outVars) {
    //potentially add Rcut
    outVars.add(ID,&Rs,1);
  }

  // right now this only does a numerical two body short range jastrow
  // based on the short range part from the breakup handled by the LRHandler
  OrbitalBase* RPAPBCConstraints::createSRTwoBody(ParticleSet& target) {
    typedef SplineJastrow<RealType> FuncType;
    typedef LinearGrid<RealType> GridType;
    
    //setRadialGrid(target);
    HandlerType* handler = LRJastrowSingleton::getHandler(target);
    RealType Rcut = 0.999999999 * handler->Basis.get_rc();
    myGrid = new GridType;
    myGrid->set(0,Rcut,51);

    /*
    for (int i = 0; i < myGrid->size(); i++) {
      RealType r=(*myGrid)(i);
      cout << r << "   " << handler->evaluate(r,1./r) << "   " << handler->evaluateLR(r) << endl;
    }
    */
  
    //create the numerical functor
    FuncType* nfunc = new FuncType;
    ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(handler);
    nfunc->initialize(sra, myGrid);
    
    TwoBodyJastrowOrbital<FuncType> *J2 = new TwoBodyJastrowOrbital<FuncType>(target);
    for (int i=0; i<4; i++) J2->addFunc(nfunc);

    return J2;
  }
     
  OrbitalBase* RPAPBCConstraints::createLRTwoBody(ParticleSet& target) {
    HandlerType* handler = LRJastrowSingleton::getHandler(target);
    LRTwoBodyJastrow* LROrbital = new LRTwoBodyJastrow(target, handler);
    return LROrbital;
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
