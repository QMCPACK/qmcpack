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

  //@{implementation of the base class of DiffOrbitalBase
  void DiffOrbitalBase::resize(int nptcls)
  {
    gradLogPsi.resize(nptcls);
    lapLogPsi.resize(nptcls);
  }

  void DiffOrbitalBase::setParameter(const string& a, RealType v)
  {
    targetParam.first=a;
    targetParam.second=v;
  }

  bool DiffOrbitalBase::isOptimizable(OptimizableSetType& optVariables)
  {
    return optVariables.find(targetParam.first) != optVariables.end();
  }

  void DiffOrbitalBase::evaluateDerivatives(ParticleSet& P, RealType ke0)
  {
    gradLogPsi=0.0;
    lapLogPsi=0.0;

    dLogPsi=differentiate(P);
    dHPsi=-0.5*Sum(lapLogPsi)-Dot(P.G,gradLogPsi)-dLogPsi*ke0;
  }
  //@}

  //@{implementation of NumericalDiffOrbital
  void NumericalDiffOrbital::resetParameters(OptimizableSetType& optVariables)
  {
    OptimizableSetType::iterator oit(optVariables.find(targetParam.first));
    if(oit!=optVariables.end())
    {
      targetParam.second=(*oit).second;
    }
  }

  void NumericalDiffOrbital::resetTargetParticleSet(ParticleSet& P)
  {
    //do nothing
  }

  DiffOrbitalBase::RealType NumericalDiffOrbital::differentiate(ParticleSet& P)
  {
    OptimizableSetType v;
    RealType curvar=targetParam.second;

    GradVectorType dg_p(gradLogPsi.size());
    GradVectorType dg_m(gradLogPsi.size());
    ValueVectorType dl_p(gradLogPsi.size());
    ValueVectorType dl_m(gradLogPsi.size());
    
    //this should be modified later
    const RealType delta=0.001;

    //reset reference orbital with +
    v[targetParam.first]=curvar+delta;
    refOrbital->resetParameters(v);
    RealType plus=refOrbital->evaluateLog(P,dg_p,dl_p);

    //reset reference orbital with -
    v[targetParam.first]=curvar-delta;
    refOrbital->resetParameters(v);
    RealType minus=refOrbital->evaluateLog(P,dg_m,dl_m);

    //restore the variable to the original state
    v[targetParam.first]=curvar;
    refOrbital->resetParameters(v);
    
    const RealType dh=1.0/(2.0*delta);
    gradLogPsi=dh*(dg_p-dg_m);
    lapLogPsi=dh*(dl_p-dl_m);
    return dh*(plus-minus);
  }
  //@}
  
  //@{implementation of AnalyticDiffOrbital
  void AnalyticDiffOrbital::resetParameters(OptimizableSetType& optVariables)
  {
    OptimizableSetType::iterator oit(optVariables.find(targetParam.first));
    if(oit!=optVariables.end())
    {
      targetParam.second=(*oit).second;
      refOrbital->resetParameters(optVariables);
    }
  }

  void AnalyticDiffOrbital::resetTargetParticleSet(ParticleSet& P)
  {
    refOrbital->resetTargetParticleSet(P);
  }

  DiffOrbitalBase::RealType AnalyticDiffOrbital::differentiate(ParticleSet& P)
  {
    return refOrbital->evaluateLog(P,gradLogPsi,lapLogPsi);
  }
  //@}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/

