//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#include "QMCDrivers/QMCUpdateBase.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommCreate.h"
#include "Message/OpenMP.h"

namespace qmcplusplus { 

  /// Constructor.
  QMCUpdateBase::QMCUpdateBase(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): W(w),Psi(psi),H(h), UpdatePbyP(true),
      RandomGen(rg), MaxAge(0),  branchEngine(0), Estimators(0)
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
        , compEstimator(0)
#endif
      { }

  /// destructor
  QMCUpdateBase::~QMCUpdateBase() 
  { 
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator) delete compEstimator;
#endif
  }

  bool QMCUpdateBase::put(xmlNodePtr cur)
  {
    //nonlocal operator is very light
    bool s= nonLocalOps.put(cur);
    s=myParams.put(cur);

#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator == 0)
    {
      //check if estimator needs to be constructed
      cur=cur->children;
      vector<string> elist;
      while(cur != NULL)
      {
        string cname((const char*)(cur->name));
        if(cname == "estimator") {
          string ename("0");
          OhmmsAttributeSet att;
          att.add(ename,"name");
          att.put(cur);
          if(ename == "gofr" || ename == "sk") //only accept gofr/sk
          {
            elist.push_back(ename);
          }
        }
        cur=cur->next;
      }
      if(elist.size())
      {
        compEstimator = new CompositeEstimatorSet(W);
      }
    }

    if(compEstimator)  compEstimator->open(-1);
#endif
    return s;
  }

  void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManager* est) {

    Estimators=est;
    branchEngine=brancher;
    branchEngine->setEstimatorManager(est);

    NumPtcl = W.getTotalNum();
    deltaR.resize(NumPtcl);
    drift.resize(NumPtcl);

    Tau=brancher->Tau;
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = std::sqrt(Tau);
  }

  void QMCUpdateBase::resetEtrial(RealType et) {
    branchEngine->E_T=et;
    branchEngine->flush(0);
  }

  void QMCUpdateBase::startRun(int blocks, bool record) 
  {
    Estimators->start(blocks,record);
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator) {
      compEstimator->open(-1);
    }
#endif
  }

  void QMCUpdateBase::stopRun() 
  {
    Estimators->stop();
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator) compEstimator->close();
#endif
  }

  void QMCUpdateBase::startBlock(int steps) {
    Estimators->startBlock(steps);
    nAccept = 0; 
    nReject=0;
    nAllRejected=0;
    nNodeCrossing=0;
    NonLocalMoveAccepted=0;
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator) compEstimator->startBlock(steps);
#endif
  }

  void QMCUpdateBase::stopBlock() {
    Estimators->stopBlock(acceptRatio());
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(compEstimator) compEstimator->stopBlock(-1,-1);
#endif
  }

  void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {
    UpdatePbyP=false;
    while(it != it_end) {
      W.R = (*it)->R;
      W.update();
      RealType logpsi(Psi.evaluateLog(W));
      setScaledDrift(Tau,W.G,(*it)->Drift);
      RealType ene = H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),ene);
      H.saveProperty((*it)->getPropertyBase());
      ++it;
    }
  }

  void QMCUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end) 
  {
    UpdatePbyP=true;
    NumPtcl=(*it)->size();//resize it always
    G.resize(NumPtcl);
    dG.resize(NumPtcl);
    L.resize(NumPtcl);
    dL.resize(NumPtcl);

    while(it != it_end) {
      (*it)->DataSet.clear();
      (*it)->DataSet.rewind();
      W.registerData(**it,(*it)->DataSet);
      RealType logpsi=Psi.registerData(W,(*it)->DataSet);

      //RealType scale=getDriftScale(Tau,W.G);
      //(*it)->Drift = scale*W.G;
      setScaledDrift(Tau,W.G,(*it)->Drift);

      RealType ene = H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),ene);
      H.saveProperty((*it)->getPropertyBase());

      ++it;
    } 
  }

  void QMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end) {
    while(it != it_end) {
      Walker_t::Buffer_t& w_buffer((*it)->DataSet);
      w_buffer.rewind();
      W.updateBuffer(**it,w_buffer);
      RealType logpsi=Psi.updateBuffer(W,w_buffer);
      RealType enew= H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),enew);
      H.saveProperty((*it)->getPropertyBase());

      //RealType scale=getDriftScale(Tau,W.G);
      //(*it)->Drift = scale*W.G;
      setScaledDrift(Tau,W.G,(*it)->Drift);

      ++it;
    }
  }

  void QMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end) {
    while(it != it_end) {
      RealType M=(*it)->Weight;
      if((*it)->Age>MaxAge) 
        M = std::min(0.5,M);
      else if((*it)->Age > 0) 
        M = std::min(1.0,M);
      (*it)->Multiplicity = M + RandomGen();
      ++it;
    }
  }

  void QMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip) {
    char fname[16];
    sprintf(fname,"test.%i",ip);
    ofstream fout(fname,ios::app);
    int i=0;
    fout << "benchMark started." << endl;
    while(it != it_end) {
      Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen); 
      W.R = m_sqrttau*deltaR+ thisWalker.R;
      W.update();
      ValueType logpsi(Psi.evaluateLog(W));
      RealType e = H.evaluate(W);
      fout << W.R[0] << W.G[0] << endl;
      fout <<  i << " " << logpsi << " " << e << endl;
      ++it;++i;
    }
    fout << "benchMark completed." << endl;
  }

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1618 $   $Date: 2007-01-14 18:10:10 -0600 (Sun, 14 Jan 2007) $
 * $Id: QMCUpdateBase.cpp 1618 2007-01-15 00:10:10Z jnkim $
 ***************************************************************************/
