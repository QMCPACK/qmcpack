//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DMCPSIVALUE_H
#define QMCPLUSPLUS_DMCPSIVALUE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

struct DMCPsiValue: public QMCHamiltonianBase
{
  typedef std::map<std::string,ParticleSet*> PtclPoolType;

  TrialWaveFunction* trialPsi;
  int psiindx,Eindx,blockSeries,counter;
  Return_t Tau,E0,EW,W,PW,Ps0;

  DMCPsiValue( )
  {
    UpdateMode.set(OPTIMIZABLE,1);
    ///defaults
    Tau=0.01;
    counter=1;
    E0=0.0;
    blockSeries=100;
  }
  ///destructor
  ~DMCPsiValue() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    trialPsi->resizeTempP(P);
  }

  void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex=plist.add("Psi");
    plist.add("PsiTW");
    plist.add("E_L_W");
    plist.add("W");
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Ps0;
    plist[myIndex+1]=PW;
    plist[myIndex+2]=EW;
    plist[myIndex+3]=W;
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    plist[myIndex+offset]=Ps0;
    plist[myIndex+1+offset]=PW;
    plist[myIndex+2+offset]=EW;
    plist[myIndex+3+offset]=W;
  }



  inline Return_t
  rejectedMove(ParticleSet& P)
  {
    evaluate(P);
    return 0.0;
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    if ( (tWalker->PropertyHistory[psiindx])[0] == 0)
    {
      Ps0=std::exp(P.Properties(LOGPSI));
      tWalker->addPropertyHistoryPoint(psiindx, Ps0 );
    }
    Return_t tEsum = std::exp(-1.0*Tau*(tWalker->getPropertyHistorySum(Eindx,blockSeries) ) );
    tWalker->addPropertyHistoryPoint(Eindx,  tWalker->Properties(LOCALENERGY)-E0 );
    Ps0=(tWalker->PropertyHistory[psiindx][0]);
    PW= tEsum *Ps0;
    EW= tEsum *( tWalker->Properties(LOCALENERGY));
    W= tEsum ;
    return 0.0;
  }

  inline Return_t
  evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t
  registerData(ParticleSet& P, BufferType& buffer)
  {
    evaluate(P);
    buffer.add(0.0);
    return 0.0;
  }

  inline Return_t
  updateBuffer(ParticleSet& P, BufferType& buffer)
  {
    evaluate(P);
    return 0.0;
  }

  inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.get(Value);
  }

  inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.put(0.0);
  }

  inline Return_t
  evaluatePbyP(ParticleSet& P, int active)
  {
    evaluate(P);
    return 0.0;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool put(xmlNodePtr cur,ParticleSet* qp,
           PtclPoolType& pset, Communicate* c )
  {
    ///get information for length of energy vector
    OhmmsAttributeSet Tattrib;
    Tattrib.add(blockSeries,"length");
    Tattrib.add(Tau,"tau");
    Tattrib.add(E0,"E0");
    Tattrib.put(cur);
    xmlNodePtr tcur = cur->children;
    ///Making wavefunction for VMC overlap calculation
    WaveFunctionFactory* WFF = new WaveFunctionFactory(qp,pset,c);
    WFF->setReportLevel(0);
    WFF->setName("DMCPsiVal");
    WFF->build(cur,false);
    trialPsi = (WFF->targetPsi)->makeClone(*qp);
    trialPsi->resizeTempP(*qp);
    delete WFF;
    ///Create a energy history vector to facilitate the DMC correction
    psiindx = qp->addPropertyHistory(1);
    Eindx = qp->addPropertyHistory(blockSeries);
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "DMCPsiValue";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    DMCPsiValue* acopy = new DMCPsiValue( );
    acopy->trialPsi = trialPsi->makeClone(qp);
    acopy->trialPsi->resizeTempP(qp);
    acopy->Eindx = Eindx;
    acopy->psiindx= psiindx;
    acopy->Tau=Tau;
    acopy->E0=E0;
    acopy->blockSeries=blockSeries;
    acopy->counter=1;
    return acopy;
  }

};
}
#endif


