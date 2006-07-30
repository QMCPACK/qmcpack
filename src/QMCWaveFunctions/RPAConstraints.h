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
#include "QMCWaveFunctions/NumericalJastrowFunctor.h"
#include "QMCWaveFunctions/RPAJastrow.h"

namespace qmcplusplus {

  struct RPAConstraints: public OrbitalConstraintsBase {
    typedef RPAJastrow<RealType> FuncType;
    bool IgnoreSpin;
    RealType Rs;
    string ID;
    vector<FuncType*> FuncList;

    ~RPAConstraints();

    RPAConstraints(bool nospin=true):IgnoreSpin(nospin) {}

    inline void apply() {
      for(int i=0; i<FuncList.size(); i++) {
        FuncList[i]->reset(Rs);
      }
    }

    void addOptimizables(VarRegistry<RealType>& outVars) {
      outVars.add(ID,&Rs,1);
    }

    OrbitalBase* createTwoBody(ParticleSet& target);
    OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
    bool put(xmlNodePtr cur);
  };

  struct RPAPBCConstraints: public OrbitalConstraintsBase {
    ///analytic functor
    typedef RPAJastrow<RealType> InFuncType;
    ///numerical functor
    typedef NumericalJastrow<RealType> FuncType;
    bool IgnoreSpin;
    RealType Rs;
    string ID;
    vector<InFuncType*> InFuncList;
    vector<FuncType*> FuncList;

    ~RPAPBCConstraints();
    RPAPBCConstraints(bool nospin=true):IgnoreSpin(nospin) {}
    void apply();
    void addOptimizables(VarRegistry<RealType>& outVars);
    OrbitalBase* createTwoBody(ParticleSet& target);
    OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
    bool put(xmlNodePtr cur);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
