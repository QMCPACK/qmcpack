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
#ifndef QMCPLUSPLUS_DMCPSIVALUE_H
#define QMCPLUSPLUS_DMCPSIVALUE_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  struct DMCPsiValue: public QMCHamiltonianBase {
  typedef map<string,ParticleSet*> PtclPoolType;
  
  TrialWaveFunction* trialPsi;
  int psiindx,Eindx,blockSeries,counter;
  Return_t Tau,E0,EW,W,PW;

    DMCPsiValue( ){
      UpdateMode.set(OPTIMIZABLE,1);
      ///defaults
       Tau=0.01;
       counter=1;
       E0=0.0;
       blockSeries=100; 
       
    }
    ///destructor
    ~DMCPsiValue() { }

    void resetTargetParticleSet(ParticleSet& P) {
      trialPsi->resizeTempP(P);
    }
    void addObservables(PropertySetType& plist)
    {
      myIndex=plist.add("E_L_W");
      plist.add("PsiTW");
      plist.add("W");
    }

    void setObservables(PropertySetType& plist)
    {
      plist[myIndex]=EW;
      plist[myIndex+1]=PW;
      plist[myIndex+2]=W;
    }
    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
      plist[myIndex+offset]=EW;
      plist[myIndex+1+offset]=PW;
      plist[myIndex+2+offset]=W;
    }
    
    
    
    inline Return_t 
    rejectedMove(ParticleSet& P) {
//         vector<double>::iterator Vit=Values.begin();
//         for(int i=0;i<nObservables;i++){
// 	  for(int j=0;j<walkerLengths[i].size();j++,Vit++){
//             (*Vit) = tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
//           }
//         }
// 	std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
      tWalker->rejectedMove();
      return 0.0;
    }
    
    inline Return_t 
    evaluate(ParticleSet& P) {
      if ( (tWalker->PropertyHistory[psiindx])[0] ==0){
        tWalker->addPropertyHistoryPoint(psiindx, std::exp( tWalker->Properties(LOGPSI)) );
      }
      
//       app_log()<<tWalker->Properties(LOCALENERGY)-E0<<"  "<<-1.0*Tau*tWalker->getPropertyHistorySum(Eindx,blockSeries)<<endl;
      Return_t tEsum = std::exp(-1.0*Tau*(tWalker->getPropertyHistorySum(Eindx,blockSeries)-0.5*tWalker->Properties(LOCALENERGY)) );
//       app_log()<<tEsum<<"  "<<tWalker->Weight<<endl;
      tWalker->addPropertyHistoryPoint(Eindx,  tWalker->Properties(LOCALENERGY)-E0 );
//       app_log()<<tEsum * ( tWalker->PropertyHistory[psiindx][0])<<endl;
//       Value= tEsum * ( tWalker->PropertyHistory[psiindx][0]);
      PW= tEsum *( tWalker->PropertyHistory[psiindx][0]);
      EW= tEsum *( tWalker->Properties(LOCALENERGY));
      W=tEsum*tWalker->Weight;
//       Value= tEsum * P.PropertyList[LOCALENERGY];
//       tWalker->Weight = tEsum;
      return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    inline Return_t 
    registerData(ParticleSet& P, BufferType& buffer) 
    {
      if ( (tWalker->PropertyHistory[psiindx])[0] ==0){
        tWalker->addPropertyHistoryPoint(psiindx, std::exp( tWalker->Properties(LOGPSI)) );
      }
      tWalker->addPropertyHistoryPoint(Eindx,  tWalker->Properties(LOCALENERGY)-E0 );
      Return_t tEsum = std::exp(-1.0*Tau*tWalker->getPropertyHistorySum(Eindx,blockSeries));
//       app_log()<<tEsum<<"  "<<tWalker->Properties(LOCALENERGY)<<endl;
//       app_log()<<tEsum * ( tWalker->PropertyHistory[psiindx][0])<<endl;
//       Value= tEsum * ( tWalker->PropertyHistory[psiindx][0]);
      Value= tEsum *( tWalker->PropertyHistory[psiindx][0]);
      
      Value= tEsum * P.PropertyList[LOCALENERGY];
      
      buffer.add(Value);
      return Value;
    }

    inline Return_t 
    updateBuffer(ParticleSet& P, BufferType& buffer) 
    {
      if ( (tWalker->PropertyHistory[psiindx])[0] ==0){
        tWalker->addPropertyHistoryPoint(psiindx, std::exp( tWalker->Properties(LOGPSI)) );
      }
      tWalker->addPropertyHistoryPoint(Eindx,  tWalker->Properties(LOCALENERGY)-E0 );
      Return_t tEsum = std::exp(-1.0*Tau*tWalker->getPropertyHistorySum(Eindx,blockSeries));
//       app_log()<<tEsum<<"  "<<tWalker->Properties(LOCALENERGY)<<endl;
//       app_log()<<tEsum * ( tWalker->PropertyHistory[psiindx][0])<<endl;
//       Value= tEsum * ( tWalker->PropertyHistory[psiindx][0]);
      Value= tEsum *( tWalker->PropertyHistory[psiindx][0]);
      
      Value= tEsum * P.PropertyList[LOCALENERGY];
      
      return Value;
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
      if ( (tWalker->PropertyHistory[psiindx])[0] ==0){
        tWalker->addPropertyHistoryPoint(psiindx, std::exp( tWalker->Properties(LOGPSI)) );
      }
      tWalker->addPropertyHistoryPoint(Eindx,  tWalker->Properties(LOCALENERGY)-E0 );
      Return_t tEsum = std::exp(-1.0*Tau*tWalker->getPropertyHistorySum(Eindx,blockSeries));
//       app_log()<<tEsum<<"  "<<tWalker->Properties(LOCALENERGY)<<endl;
//       app_log()<<tEsum * ( tWalker->PropertyHistory[psiindx][0])<<endl;
//       Value= tEsum * ( tWalker->PropertyHistory[psiindx][0]);
      Value= tEsum *( tWalker->PropertyHistory[psiindx][0]);
      
      Value= tEsum * P.PropertyList[LOCALENERGY];
      
      return Value;
    }
 
 bool put(xmlNodePtr cur) {
   return true;
    }
    
     bool put(xmlNodePtr cur,ParticleSet* qp, 
    PtclPoolType& pset, Communicate* c ) {
       ///get information for length of energy vector

       OhmmsAttributeSet Tattrib; 
       Tattrib.add(blockSeries,"length");
       Tattrib.add(Tau,"tau");
       Tattrib.add(E0,"E0");
       Tattrib.put(cur);
       
       xmlNodePtr tcur = cur->children;
//       app_log()<<"PUTTING LOGPSI "<<((const char*)(tcur->name))<<endl;
      
      ///Making wavefunction for VMC overlap calculation
      WaveFunctionFactory* WFF = new WaveFunctionFactory(qp,pset,c);
      WFF->setReportLevel(0);
      WFF->setName("HamPsiVal");
      WFF->build(cur,false);
      
      trialPsi = (WFF->targetPsi)->makeClone(*qp);
      
      trialPsi->resizeTempP(*qp);
      
      delete WFF;
      
      ///Create a energy history vector to facilitate the DMC correction
      Eindx = qp->addPropertyHistory(blockSeries);
      psiindx = qp->addPropertyHistory(1);
      
      return true;
    }

    bool get(std::ostream& os) const {
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

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3111 $   $Date: 2008-09-14 13:20:23 -0500 (Sun, 14 Sep 2008) $
 * $Id: PsiValue.h 3111 2008-09-14 18:20:23Z jnkim $ 
 ***************************************************************************/

