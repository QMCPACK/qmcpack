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
#include "Message/OpenMP.h"

namespace qmcplusplus
  {

  /// Constructor.
  QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, RandomGenerator_t& rg)
      : W(w),Psi(psi),Guide(guide),H(h)
      , UpdatePbyP(true), UseTMove(false)
      , NumPtcl(0), nSubSteps(1)
      , RandomGen(rg), MaxAge(0),  m_r2max(-1)
      , branchEngine(0), Estimators(0)
  {
    myParams.add(m_r2max,"maxDisplSq","double"); //maximum displacement
    myParams.add(nSubSteps,"subSteps","int");
    myParams.add(nSubSteps,"substeps","int");
    myParams.add(nSubSteps,"sub_steps","int");
  }

  /// Constructor.
  QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
      : W(w),Psi(psi),H(h),Guide(psi)
      , UpdatePbyP(true), UseTMove(false)
      , NumPtcl(0), nSubSteps(1)
      , RandomGen(rg), MaxAge(0),  m_r2max(-1)
      , branchEngine(0), Estimators(0)
  {
    myParams.add(m_r2max,"maxDisplSq","double"); //maximum displacement
    myParams.add(nSubSteps,"subSteps","int");
    myParams.add(nSubSteps,"substeps","int");
    myParams.add(nSubSteps,"sub_steps","int");
  }
  
  /// destructor
  QMCUpdateBase::~QMCUpdateBase()
  {
  }

  bool QMCUpdateBase::put(xmlNodePtr cur)
  {
    //nonlocal operator is very light
    UseTMove = nonLocalOps.put(cur);
    bool s=myParams.put(cur);
    if (branchEngine) branchEngine->put(cur);
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

    G.resize(NumPtcl);
    dG.resize(NumPtcl);
    L.resize(NumPtcl);
    dL.resize(NumPtcl);

    ////Tau=brancher->getTau();
    ////m_oneover2tau = 0.5/Tau;
    ////m_sqrttau = std::sqrt(Tau);
    SpeciesSet tspecies(W.getSpeciesSet());
    int massind=tspecies.addAttribute("mass");
    RealType mass = tspecies(massind,0);
    //if(mass<numeric_limits<RealType>::epsilon())
    //{
    //  mass=1.0;
    //  tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]))=1.0;
    //}

    //oneovermass is NOT a data member and used only here!!!
    //use m_ if using all lowercase for the variables
    RealType oneovermass = 1.0/mass;
    RealType oneoversqrtmass = std::sqrt(oneovermass);
    Tau=brancher->getTau();
    m_tauovermass = Tau/mass;
    m_oneover2tau = 0.5/(m_tauovermass);
    m_sqrttau = std::sqrt(m_tauovermass);

    //app_log() << "  QMCUpdateBase::resetRun m/tau=" << m_tauovermass << endl;
    if (m_r2max<0)
      m_r2max =  W.Lattice.LR_rc* W.Lattice.LR_rc;

    //app_log() << "  Setting the bound for the displacement max(r^2) = " <<  m_r2max << endl;
  }

  void QMCUpdateBase::resetEtrial(RealType et)
  {
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

  void QMCUpdateBase::startBlock(int steps)
  {
    Estimators->startBlock(steps);
    nAccept = 0;
    nReject=0;
    nAllRejected=0;
    nNodeCrossing=0;
    NonLocalMoveAccepted=0;
  }

  void QMCUpdateBase::stopBlock(bool collectall)
  {
    Estimators->stopBlock(acceptRatio(),collectall);
  }

  void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
  {
    UpdatePbyP=false;
    for (;it != it_end; ++it)
      {
        W.R = (*it)->R;
        W.update();
        RealType logpsi(Psi.evaluateLog(W));
        //setScaledDriftPbyP(Tau*m_oneovermass,W.G,(*it)->Drift);
        RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
        RealType ene = H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0, nodecorr);
        (*it)->Weight=1;
        H.saveProperty((*it)->getPropertyBase());
      }
  }

  void QMCUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
  {
    UpdatePbyP=true;

    for (;it != it_end; ++it)
      {
        Walker_t& thisWalker(**it);
        W.loadWalker(thisWalker,UpdatePbyP);
        if (thisWalker.DataSet.size())
          thisWalker.DataSet.clear();
        
        Walker_t::Buffer_t tbuffer;
        RealType logpsi=Psi.registerData(W,tbuffer);
        thisWalker.DataSet=tbuffer;

        //setScaledDriftPbyP(m_tauovermass,W.G,(*it)->Drift);
        RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
        RealType ene = H.evaluate(W);

        thisWalker.resetProperty(logpsi,Psi.getPhase(),ene, 0.0,0.0, nodecorr);
        H.saveProperty(thisWalker.getPropertyBase());
        thisWalker.ReleasedNodeAge=0;
        thisWalker.ReleasedNodeWeight=0;
        thisWalker.Weight=1;
      }
  }

  QMCUpdateBase::RealType
  QMCUpdateBase::getNodeCorrection(const ParticleSet::ParticleGradient_t& g, ParticleSet::ParticlePos_t& gscaled)
  {
//       PAOps<RealType,OHMMS_DIM>::copy(g,gscaled);
    //// DriftOperators.h getNodeCorrectionP
    //RealType norm=0.0, norm_scaled=0.0;
    //for(int i=0; i<g.size(); ++i)
    //{
    //  RealType vsq=dot(g[i],g[i]);
    //  RealType x=vsq*Tau;
    //  RealType scale= (vsq<numeric_limits<RealType>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
    //  norm_scaled+=vsq*scale*scale;
    //  norm+=vsq;
    //}
    //return std::sqrt(norm_scaled/norm);

    // DriftOperators.h getNodeCorrectionW
    setScaledDrift(m_tauovermass,g,gscaled);
    RealType vsq=Dot(g,g);
    RealType x=m_tauovermass*vsq;
    return (vsq<numeric_limits<RealType>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
  }

  void QMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end)
  {

    for (;it != it_end; ++it)
      {
        Walker_t& thisWalker(**it);
        W.loadWalker(thisWalker,UpdatePbyP);
        Walker_t::Buffer_t& w_buffer((*it)->DataSet);

        RealType logpsi=Psi.updateBuffer(W,w_buffer,true);

        //needed to copy R/L/G
        W.saveWalker(thisWalker);
//         thisWalker.Weight=1;
        //thisWalker.Properties(DRIFTSCALE)=getNodeCorrection(W.G,(*it)->Drift);
        //RealType enew= H.evaluate(W);
        //thisWalker.resetProperty(logpsi,Psi.getPhase(),enew);
        //H.saveProperty(thisWalker.getPropertyBase());

        //(*it)->Drift=W.G;//copy gradients to drift
        //scaling factor per particle
        //setScaledDriftPbyP(Tau,W.G,(*it)->Drift);

        ////calculate the scaling factor
        //RealType scale=getDriftScale(Tau,W.G);
        //assignDrift(scale,W.G,(*it)->Drift);

        ////This is the original
        //setScaledDrift(Tau,W.G,(*it)->Drift);
      }
  }

  void QMCUpdateBase::setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
  {
    for (;it != it_end; ++it)
      {
        RealType M=std::abs((*it)->Weight);
        (*it)->Multiplicity = std::floor(M + RandomGen());
      }
  }

  void QMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
  {
    for (;it != it_end; ++it)
      {
        RealType M=(*it)->Weight;
        if ((*it)->Age>MaxAge)
          M = std::min(0.5,M);
        else if ((*it)->Age > 0)
          M = std::min(1.0,M);
        (*it)->Multiplicity = M + RandomGen();
      }
  }

  void QMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip)
  {
    char fname[16];
    sprintf(fname,"test.%i",ip);
    ofstream fout(fname,ios::app);
    int i=0;
    fout << "benchMark started." << endl;
    for (;it != it_end; ++it,++i)
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
