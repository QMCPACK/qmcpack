//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#define QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus {

  class LRTwoBodyJastrow;

  struct RPAPBCConstraints: public OrbitalConstraintsBase {

    typedef LRHandlerBase HandlerType;    
    
    enum {USE_BREAKUP=0, USE_RPA, USE_NOME};

    bool IgnoreSpin;
    bool DropLongRange;
    bool DropShortRange;
    int LongRangeForm;
    RealType Rs;
    RealType Kc;
    string ID_Rs;
    string MyName;

    RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

    ~RPAPBCConstraints();
    void addOptimizables(OptimizableSetType& outVars);
    void resetParameters(OptimizableSetType& optVariables);
    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);
    void addExtra2ComboOrbital(ComboOrbital* jcombo);
    bool put(xmlNodePtr cur);

  private:
    /** main handler
     */
    HandlerType* myHandler;
    /** handler which performs a breakup in real and k-space
     *
     * ShortRangeRPA uses realHanlder
     */
    HandlerType* realHandler;
    /** hanlder which can overwrite the long-range part of a RPA jastro
     */
    HandlerType* kspaceHandler;

    LRTwoBodyJastrow* LongRangeRPA;
    OrbitalBase* ShortRangeRPA;
  };

 // struct RPAConstraints: public OrbitalConstraintsBase {
 //   typedef RPAJastrow<RealType> FuncType;
 //   bool IgnoreSpin;
 //   RealType Rs;
 //   string ID_Rs;
 //   vector<FuncType*> FuncList;

 //   ~RPAConstraints();

 //   RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

 //   void addOptimizables(OptimizableSetType& outVars); 
 //   /** update the optimizable variables
 //    */
 //   void resetParameters(OptimizableSetType& optVariables);


 //   OrbitalBase* createTwoBody();
 //   OrbitalBase* createOneBody(ParticleSet& source);
 //   void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
 //   bool put(xmlNodePtr cur);
 // };
     
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
