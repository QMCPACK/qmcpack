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
#ifndef QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#define QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/ComboOrbital.h"

namespace qmcplusplus {

  struct WMConstraints: public OrbitalConstraintsBase {
    ///basis type
    typedef OptimizableFunctorBase<RealType>  BasisType;
    ///basis set type
    typedef vector<BasisType*> BasisSetType;
    ///analytic functor
    typedef ComboFunctor<RealType> InFuncType;
    ///numerical functor
    typedef CubicBsplineSingle<RealType> FuncType;
    ///flag to tunr on/off spin-dependent term, always off
    bool IgnoreSpin;
    ///cutoff radius
    RealType Rcut;

    map<string,BasisSetType*> myBasisSet;
    vector<InFuncType*> InFuncList;
    vector<FuncType*> FuncList;

    ~WMConstraints();
    WMConstraints(bool nospin=true):IgnoreSpin(nospin) {}
    void apply();
    void addOptimizables(VarRegistry<RealType>& outVars);
    OrbitalBase* createTwoBody(ParticleSet& target);
    OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
    inline void addTwoBodyPart(ParticleSet& target, ComboOrbital* jcombo) {
      OrbitalBase* twobody = createTwoBody(target);
      if (twobody) jcombo->Psi.push_back(twobody);
    }
    bool put(xmlNodePtr cur);
    void addBasisGroup(xmlNodePtr cur);
    InFuncType* createCorrelation(xmlNodePtr cur, BasisSetType* basis);
    InFuncType* createDefaultTwoBody(xmlNodePtr cur, const string& tname);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
