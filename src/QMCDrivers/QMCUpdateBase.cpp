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
  QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg): 
    W(w),Psi(psi),H(h), UpdatePbyP(true),
    RandomGen(rg), MaxAge(0),  m_r2max(-1), branchEngine(0), Estimators(0)
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
      , compEstimator(0)
#endif
  { 
    myParams.add(m_r2max,"maxDisplSq","double"); //maximum displacement
  }

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

  void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManager* est) 
  {

    Estimators=est;
    branchEngine=brancher;
    branchEngine->setEstimatorManager(est);

    NumPtcl = W.getTotalNum();
    deltaR.resize(NumPtcl);
    drift.resize(NumPtcl);

    Tau=brancher->getTau();
    m_oneover2tau = 0.5/Tau;
    m_sqrttau = std::sqrt(Tau);

    if(m_r2max<0)
      m_r2max = W.Lattice.LR_rc* W.Lattice.LR_rc;

    //app_log() << "  Setting the bound for the displacement max(r^2) = " <<  m_r2max << endl;
  }

  void QMCUpdateBase::resetEtrial(RealType et) {
    //branchEngine->E_T=et;
    branchEngine->setTrialEnergy(et,1.0);
    branchEngine->flush(0);
  }

  void QMCUpdateBase::startRun(int blocks, bool record) 
  {
    Estimators->start(blocks,record);
  }

  void QMCUpdateBase::stopRun() 
  {
    Estimators->stop();
  }

  void QMCUpdateBase::startBlock(int steps) {
    Estimators->startBlock(steps);
    nAccept = 0; 
    nReject=0;
    nAllRejected=0;
    nNodeCrossing=0;
    NonLocalMoveAccepted=0;
  }

  void QMCUpdateBase::stopBlock() {
    Estimators->stopBlock(acceptRatio());
  }

  void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {
    UpdatePbyP=false;
    for(;it != it_end; ++it)
    {
      W.R = (*it)->R;
      W.update();
      RealType logpsi(Psi.evaluateLog(W));
      setScaledDrift(Tau,W.G,(*it)->Drift);
      RealType ene = H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),ene);
      H.saveProperty((*it)->getPropertyBase());
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

    for(;it != it_end; ++it)
    {
      Walker_t::Buffer_t tbuffer;
      W.registerData(**it,tbuffer);
      RealType logpsi=Psi.registerData(W,tbuffer);
      (*it)->DataSet=tbuffer;

      //RealType scale=getDriftScale(Tau,W.G);
      //(*it)->Drift = scale*W.G;
      setScaledDrift(Tau,W.G,(*it)->Drift);

      RealType ene = H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),ene);
      H.saveProperty((*it)->getPropertyBase());
    } 

  }

  void QMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end) {

    RealType tauinv=1.0/Tau;
    for(;it != it_end; ++it)
    {
      Walker_t::Buffer_t& w_buffer((*it)->DataSet);
      w_buffer.rewind();
      W.updateBuffer(**it,w_buffer);
      RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
      RealType enew= H.evaluate(W);

      (*it)->resetProperty(logpsi,Psi.getPhase(),enew,1.0,1.0,1.0);
      H.saveProperty((*it)->getPropertyBase());

      (*it)->Drift=W.G;//copy gradients to drift
      //scaling factor per particle
      //setScaledDriftPbyP(Tau,W.G,(*it)->Drift);
      
      ////calculate the scaling factor
      //RealType scale=getDriftScale(Tau,W.G);
      //assignDrift(scale,W.G,(*it)->Drift);
      
      ////This is the original
      //setScaledDrift(Tau,W.G,(*it)->Drift);
    }
  }

  void QMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end) {
    for(;it != it_end; ++it)
    {
      RealType M=(*it)->Weight;
      if((*it)->Age>MaxAge) 
        M = std::min(0.5,M);
      else if((*it)->Age > 0) 
        M = std::min(1.0,M);
      (*it)->Multiplicity = M + RandomGen();
    }
  }

  void QMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip) {
    char fname[16];
    sprintf(fname,"test.%i",ip);
    ofstream fout(fname,ios::app);
    int i=0;
    fout << "benchMark started." << endl;
    for(;it != it_end; ++it,++i)
    {
      Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen); 
      W.R = m_sqrttau*deltaR+ thisWalker.R;
      W.update();
      ValueType logpsi(Psi.evaluateLog(W));
      RealType e = H.evaluate(W);
      fout << W.R[0] << W.G[0] << endl;
      fout <<  i << " " << logpsi << " " << e << endl;
    }
    fout << "benchMark completed." << endl;
  }

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1618 $   $Date: 2007-01-14 18:10:10 -0600 (Sun, 14 Jan 2007) $
 * $Id: QMCUpdateBase.cpp 1618 2007-01-15 00:10:10Z jnkim $
 ***************************************************************************/
