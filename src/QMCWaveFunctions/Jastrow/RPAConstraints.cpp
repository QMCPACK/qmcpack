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
#include "QMCWaveFunctions/Jastrow/RPAConstraints.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "Utilities/IteratorUtility.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {

  ////////////////////////////////////////
  //RPAConstraints definitions
  ////////////////////////////////////////
  RPAConstraints::RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin)
    {
      JComponent.set(TWOBODY);
    }

  RPAConstraints::~RPAConstraints() {
    delete_iter(FuncList.begin(), FuncList.end());
  }

  bool RPAConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID_Rs="rs"; Rs=-1.0;
      }
    }
    else 
    {
      ID_Rs=(*vit).second.first; 
      Rs=(*vit).second.second;
    }
    return true;
  }

  void RPAConstraints::addOptimizables(OptimizableSetType& outVars) 
  {
    outVars[ID_Rs]=Rs;
  }

  void RPAConstraints::resetParameters(OptimizableSetType& optVariables) 
  {
    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
    if(it != optVariables.end()) 
    {
      Rs=(*it).second;
      for(int i=0; i<FuncList.size(); i++) {
        FuncList[i]->reset(Rs);
      }
    }
  }

  OrbitalBase* RPAConstraints::createTwoBody() {

    if(Rs<0) {
      if(targetPtcl.Lattice.SuperCellEnum) {
        Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
      } else { 
        Rs=1.0;
      }
    }

    typedef TwoBodyJastrowOrbital<FuncType> JeeType;
    JeeType *J2 = new JeeType(targetPtcl);
//    if(IgnoreSpin) {
      FuncType *func=new FuncType(false);
      func->reset(Rs);
      for(int i=0; i<4; i++) {
        J2->addFunc(func);
      }
      FuncList.push_back(func);
      app_log() << "  RPAConstraints::Adding Spin-independent RPA Two-Body Jastrow " << endl;
      app_log() << "    Rs=" << Rs << "  F=" << 1.0/func->B << endl;
//    } else {
//      FuncType *funcUU=new FuncType(true);
//      FuncType *funcUD=new FuncType(false);
//      funcUU->reset(Rs);
//      funcUD->reset(Rs);
//
//      J2->addFunc(funcUU);
//      J2->addFunc(funcUD);
//      J2->addFunc(funcUD);
//      J2->addFunc(funcUU);
//
//      FuncList.push_back(funcUU);
//      FuncList.push_back(funcUD);
//      app_log() << "  RPAConstraints::Adding Spin-dependent RPA Two-Body Jastrow " << endl;
//      app_log() << "    Rs=" << Rs << "  F(uu)=" << 1.0/funcUU->B << "  F(ud)= " << 1.0/funcUD->B << endl;
//    }
    return J2;
  }

  OrbitalBase* RPAConstraints::createOneBody(ParticleSet& source) {
    return 0;
  }

  ////////////////////////////////////////
  //RPAPBCConstraints definitions
  ////////////////////////////////////////
  RPAPBCConstraints::RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin):
    OrbitalConstraintsBase(p,psi),IgnoreSpin(nospin),LongRangeRPA(0)
    {
      JComponent.set(MULTIPLE);
      JComponent.set(TWOBODY);
      JComponent.set(LONGRANGE);
    }

  RPAPBCConstraints::~RPAPBCConstraints() 
  {  
    //may delete LongRangeRPA
  }

  bool RPAPBCConstraints::put(xmlNodePtr cur) {
    bool success=getVariables(cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find("A"));
    if(vit == inVars.end()) {
      vit=inVars.find("rs");
      if(vit == inVars.end()) {
        ID_Rs="rs"; Rs=-1.0;
      }
    }
    else
    {
      ID_Rs=(*vit).second.first; 
      Rs=(*vit).second.second;
    }
    return true;
  }

  void RPAPBCConstraints::addOptimizables(OptimizableSetType& outVars) {
    //potentially add Rcut
    outVars[ID_Rs]=Rs;
    if(LongRangeRPA) LongRangeRPA->put(NULL,outVars);
  }

  void RPAPBCConstraints::resetParameters(OptimizableSetType& optVariables) 
  { 
    OptimizableSetType::iterator it(optVariables.find(ID_Rs));
    if(it != optVariables.end()) Rs=(*it).second;
  }

  // right now this only does a numerical two body short range jastrow
  // based on the short range part from the breakup handled by the LRHandler
  OrbitalBase* RPAPBCConstraints::createSRTwoBody() {

    typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
    typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
    typedef LinearGrid<RealType> GridType;
    
    //setRadialGrid(target);
    HandlerType* handler = LRJastrowSingleton::getHandler(targetPtcl);
    //RealType Rcut = 0.999999999 * handler->Basis.get_rc();
    RealType Rcut = handler->Basis.get_rc();
    myGrid = new GridType;
    int npts=static_cast<int>(Rcut/0.05)+1;
    myGrid->set(0,Rcut,npts);
  
    //create the numerical functor
    FuncType* nfunc = new FuncType;
    ShortRangePartAdapter<RealType>* sra = new ShortRangePartAdapter<RealType>(handler);
    nfunc->initialize(sra, myGrid);

    ofstream fout("rpa.short.dat");
    for (int i = 0; i < myGrid->size(); i++) {
      RealType r=(*myGrid)(i);
      fout << r << "   " << nfunc->evaluate(r) << "   " << handler->evaluate(r,1.0/r) << " " << handler->evaluateLR(r) << endl;
    }
    
    TwoBodyJastrowOrbital<FuncType> *J2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl);
    for (int i=0; i<4; i++) J2->addFunc(nfunc);

    return J2;
  }
     
  OrbitalBase* RPAPBCConstraints::createLRTwoBody() {
    if(LongRangeRPA==0)
    {
      HandlerType* handler = LRJastrowSingleton::getHandler(targetPtcl);
      LongRangeRPA = new LRTwoBodyJastrow(targetPtcl, handler);
    }
    return LongRangeRPA;
  }

  OrbitalBase* RPAPBCConstraints::createTwoBody() 
  {
    if(Rs<0) {
      if(targetPtcl.Lattice.SuperCellEnum) {
        Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
      } else {
        Rs=1.0;
      }
    }
    app_log() << "  RPAPBCConstraints::addTwoBodyPart Rs " << Rs << endl;
    OrbitalBase* srp=createSRTwoBody();
    return srp;
  }

  void RPAPBCConstraints::addExtra2ComboOrbital(ComboOrbital* jcombo)
  {
    jcombo->Psi.push_back(createLRTwoBody());
  }

  //void RPAPBCConstraints::addTwoBodyPart(ComboOrbital* jcombo) {
  //  
  //  if(Rs<0) {
  //    if(targetPtcl.Lattice.BoxBConds[0]) {
  //      Rs=std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/static_cast<RealType>(targetPtcl.getTotalNum()),1.0/3.0);
  //    } else {
  //      Rs=1.0;
  //    }
  //  }

  //  app_log() << "  RPAPBCConstraints::addTwoBodyPart Rs " << Rs << endl;

  //  OrbitalBase* sr = createSRTwoBody();
  //  if (sr) jcombo->Psi.push_back(sr);
  //  //OrbitalBase* lr = createLRTwoBody(target);
  //  //if (lr) jcombo->Psi.push_back(lr);
  //} 
     
  OrbitalBase* RPAPBCConstraints::createOneBody(ParticleSet& source) {
    return 0;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
