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
#ifndef QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#define QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {

  class LRTwoBodyJastrow;

  struct RPAConstraints: public OrbitalConstraintsBase {
    typedef RPAJastrow<RealType> FuncType;
    bool IgnoreSpin;
    RealType Rs;
    string ID_Rs;
    vector<FuncType*> FuncList;

    ~RPAConstraints();

    RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    void addOptimizables(OptimizableSetType& outVars); 
    /** update the optimizable variables
     */
    void resetParameters(OptimizableSetType& optVariables);


    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
    bool put(xmlNodePtr cur);
  };
     
  struct RPAPBCConstraints: public OrbitalConstraintsBase {
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    
    bool IgnoreSpin;
    bool DropLongRange;
    bool DropShortRange;
    bool DummyData;
    RealType Rs;
    RealType Kc;
    string ID_Rs;

    RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    ~RPAPBCConstraints();
    void addOptimizables(OptimizableSetType& outVars);
    void resetParameters(OptimizableSetType& optVariables);
    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    void addExtra2ComboOrbital(ComboOrbital* jcombo);

    bool put(xmlNodePtr cur);
    //OrbitalBase* createSRTwoBody();
    //OrbitalBase* createLRTwoBody();
     
  private:
    HandlerType* myHandler;
    LRTwoBodyJastrow* LongRangeRPA;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
