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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/DMC/DMCUpdateBase.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/CommCreate.h"

namespace qmcplusplus {

  /// Constructor.
  DMCUpdateBase::DMCUpdateBase(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      RandomGenerator_t& rg): W(w),Psi(psi),H(h), 
                              RandomGen(rg), MaxAge(0)
    { }
  
  /// destructor
  DMCUpdateBase::~DMCUpdateBase() { }

  void DMCUpdateBase::resetRun(BranchEngineType* brancher) {
    branchEngine=brancher;
    NumPtcl = W.getTotalNum();
    deltaR.resize(NumPtcl);
    drift.resize(NumPtcl);
    
    SpeciesSet tspecies(W.getSpeciesSet());
    RealType mass = tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]));
    oneovermass = 1.0/mass;
    RealType oneoversqrtmass = std::sqrt(oneovermass);

    Tau=brancher->Tau;
    m_oneover2tau = 0.5*mass/Tau;
    m_sqrttau = sqrt(Tau*oneovermass);
  }

  void DMCUpdateBase::resetEtrial(RealType et) {
    branchEngine->E_T=et;
    branchEngine->flush(0);
  }

  void DMCUpdateBase::startBlock() {
    nAccept = 0; 
    nReject=0;
    nAllRejected=0;
    nNodeCrossing=0;
    NonLocalMoveAccepted=0;
  }
  void DMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end) {
    if(G.size() != NumPtcl) {
      G.resize(NumPtcl);
      dG.resize(NumPtcl);
      L.resize(NumPtcl);
      dL.resize(NumPtcl);
    }
    while(it != it_end) {
      (*it)->DataSet.rewind();
      W.registerData(**it,(*it)->DataSet);
      RealType logpsi=Psi.registerData(W,(*it)->DataSet);

      //RealType scale=getDriftScale(Tau,W.G);
      //(*it)->Drift = scale*W.G;
      setScaledDriftPbyP(Tau*oneovermass,W.G,(*it)->Drift);

      RealType ene = H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),ene);
      H.saveProperty((*it)->getPropertyBase());
      ++it;
    } 
  }

  void DMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end) {
    while(it != it_end) {
      Walker_t::Buffer_t& w_buffer((*it)->DataSet);
      w_buffer.rewind();
      W.updateBuffer(**it,w_buffer);
      RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
      RealType enew= H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getPhase(),enew);
      H.saveProperty((*it)->getPropertyBase());

      //RealType scale=getDriftScale(Tau,W.G);
      //(*it)->Drift = scale*W.G;
      setScaledDriftPbyP(Tau*oneovermass,W.G,(*it)->Drift);

      ++it;
    }
  }

  void DMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end) {
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

  void DMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip) {
    char fname[16];
    sprintf(fname,"test.%i",ip);
    ofstream fout(fname,ios::app);
    int i=0;
    while(it != it_end) {
      Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen); 
      W.R = m_sqrttau*deltaR+ thisWalker.R;
      W.update();
      ValueType logpsi(Psi.evaluateLog(W));
      RealType e = H.evaluate(W);
      fout <<  i << " " << logpsi << " " << e << endl;
      ++it;++i;
    }
  }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
