//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_PSIVALUE_H
#define QMCPLUSPLUS_PSIVALUE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

struct PsiValue: public QMCHamiltonianBase
{
  typedef map<string,ParticleSet*> PtclPoolType;

  TrialWaveFunction* trialPsi;
  int Power;

  PsiValue(int pwr )
  {
    Power=pwr;
    UpdateMode.set(OPTIMIZABLE,1);
  }
  ///destructor
  ~PsiValue() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    if (Power!=0)
      trialPsi->resizeTempP(P);
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    if (Power!=0)
    {
      RealType logValue = trialPsi->evaluateLogOnly(P);
      //       app_log()<<tWalker->Properties(LOGPSI)<<"  "<<logValue<<endl;
      Value = std::exp(Power*logValue);
      return Value;
    }
    else
      return Value = std::exp(tWalker->Properties(LOGPSI));
  }

  inline Return_t
  evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t
  registerData(ParticleSet& P, BufferType& buffer)
  {
    if (Power!=0)
    {
      RealType logValue = trialPsi->evaluateLogOnly(P);
      Value = std::exp(Power*logValue);
      buffer.add(Value);
      return Value;
    }
    else
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      buffer.add(Value);
      return Value;
    }
  }

  inline Return_t
  updateBuffer(ParticleSet& P, BufferType& buffer)
  {
    if (Power!=0)
    {
      RealType logValue = trialPsi->evaluateLogOnly(P);
      Value = std::exp(Power*logValue);
      buffer.add(Value);
      return Value;
    }
    else
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      buffer.add(Value);
      return Value;
    }
  }

  inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.get(Value);
  }

  inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.put(Value);
  }

  inline Return_t
  evaluatePbyP(ParticleSet& P, int active)
  {
    if (Power!=0)
    {
      RealType logValue = trialPsi->evaluateLogOnly(P);
      Value = std::exp(Power*logValue);
      return Value;
    }
    else
    {
      Value = std::exp(tWalker->Properties(LOGPSI));
      return Value;
    }
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool put(xmlNodePtr cur,ParticleSet* qp,
           PtclPoolType& pset, Communicate* c )
  {
    app_log()<<"Psi evaluator is being added with Power="<<Power<<endl;
    if (Power!=0)
    {
      xmlNodePtr tcur = cur->children;
      //       app_log()<<"PUTTING LOGPSI "<<((const char*)(tcur->name))<<endl;
      WaveFunctionFactory* WFF = new WaveFunctionFactory(qp,pset,c);
      WFF->setReportLevel(0);
      WFF->setName("HamPsiVal");
      WFF->build(cur,false);
      trialPsi = (WFF->targetPsi)->makeClone(*qp);
      trialPsi->resizeTempP(*qp);
      delete WFF;
    }
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "PsiValue";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    PsiValue* acopy = new PsiValue(Power);
    if (Power!=0)
    {
      acopy->trialPsi = trialPsi->makeClone(qp);
      acopy->trialPsi->resizeTempP(qp);
    }
    return acopy;
  }

};
}
#endif

/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 3111 $   $Date: 2008-09-14 13:20:23 -0500 (Sun, 14 Sep 2008) $
* $Id: PsiValue.h 3111 2008-09-14 18:20:23Z jnkim $
***************************************************************************/

