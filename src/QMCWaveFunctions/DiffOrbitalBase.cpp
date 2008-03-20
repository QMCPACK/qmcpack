//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "ParticleBase/ParticleAttribOps.h"

/**@file DiffOrbitalBase.cpp
 *@brief Definition of NumericalDiffOrbital
 */
namespace qmcplusplus {
  DiffOrbitalBase::DiffOrbitalBase(OrbitalBase* orb):
    FirstIndex(0),LastIndex(0)
  { 
    if(orb) refOrbital.push_back(orb);
  }

  //@{implementation of NumericalDiffOrbital
  void NumericalDiffOrbital::resetParameters(OptimizableSetType& optVars)
  {
    //do nothing
  }

  void NumericalDiffOrbital::resetTargetParticleSet(ParticleSet& P)
  {
    int nptcls=P.getTotalNum();
    dg_p.resize(nptcls);
    dl_p.resize(nptcls);
    dg_m.resize(nptcls);
    dl_m.resize(nptcls);
    gradLogPsi.resize(nptcls);
    lapLogPsi.resize(nptcls);
  }

  void NumericalDiffOrbital::evaluateDerivatives(ParticleSet& P, RealType ke0, OptimizableSetType& optVars)
  {
    OptimizableSetType v;
    RealType curvar=optVars.getValue(FirstIndex);
    string vname=optVars.getName(FirstIndex);

    //this should be modified later
    const RealType delta=0.001;

    RealType plus=0.0;
    RealType minus=0.0;
    dg_p=0.0;
    dl_p=0.0;
    dg_m=0.0;
    dl_m=0.0;

    //accumulate plus and minus displacement
    for(int i=0; i<refOrbital.size(); ++i)
    {
      //reset reference orbital with +
      v[vname]=curvar+delta;
      refOrbital[i]->resetParameters(v);
      plus+=refOrbital[i]->evaluateLog(P,dg_p,dl_p);

      //reset reference orbital with -
      v[vname]=curvar-delta;
      refOrbital[i]->resetParameters(v);
      minus+=refOrbital[i]->evaluateLog(P,dg_m,dl_m);

      //restore the variable to the original state
      v[vname]=curvar;
      refOrbital[i]->resetParameters(v);
    }
    
    const RealType dh=1.0/(2.0*delta);
    gradLogPsi=dh*(dg_p-dg_m);
    lapLogPsi=dh*(dl_p-dl_m);

    RealType dLogPsi=dh*(plus-minus);
    //optVars.setDeriv(FirstIndex,dLogPsi,-0.5*Sum(lapLogPsi)-Dot(P.G,gradLogPsi)-dLogPsi*ke0);
    optVars.setDeriv(FirstIndex,dLogPsi,-0.5*Sum(lapLogPsi)-Dot(P.G,gradLogPsi));
  }
  //@}
  
  //@{implementation of AnalyticDiffOrbital
  void AnalyticDiffOrbital::resetParameters(OptimizableSetType& optVars)
  {
    for(int i=0; i<refOrbital.size(); ++i)
      refOrbital[i]->resetParameters(optVars);
  }

  void AnalyticDiffOrbital::resetTargetParticleSet(ParticleSet& P)
  {
    for(int i=0; i<refOrbital.size(); ++i) refOrbital[i]->resetTargetParticleSet(P);
    int nptcls=P.getTotalNum();
    if(gradLogPsi.size()!=nptcls)
    {
      gradLogPsi.resize(nptcls);
      lapLogPsi.resize(nptcls);
    }
  }

  void AnalyticDiffOrbital::evaluateDerivatives(ParticleSet& P, RealType ke0, OptimizableSetType& optVars)
  {
    RealType dLogPsi=0.0;
    gradLogPsi=0.0;
    lapLogPsi=0.0;
    for(int i=0; i<refOrbital.size(); ++i)
      dLogPsi+=refOrbital[i]->evaluateLog(P,gradLogPsi,lapLogPsi);

    //optVars.setDeriv(FirstIndex,dLogPsi,-0.5*Sum(lapLogPsi)-Dot(P.G,gradLogPsi)-dLogPsi*ke0);
    optVars.setDeriv(FirstIndex,dLogPsi,-0.5*Sum(lapLogPsi)-Dot(P.G,gradLogPsi));
  }
  //@}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/

