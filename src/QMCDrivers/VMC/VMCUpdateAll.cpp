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
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus { 

  VMCUpdateAll::VMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg), nSubSteps(1)
    { 
      myParams.add(nSubSteps,"subSteps","int"); 
      myParams.add(nSubSteps,"substeps","int");
    }

  VMCUpdateAll::~VMCUpdateAll()
  {
  }

  void VMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure) 
  {
    for(; it!= it_end; ++it)
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen);
      if(!W.makeMove(thisWalker,deltaR,m_sqrttau)) continue;
      //W.R = m_sqrttau*deltaR + thisWalker.R;
      //W.update();

      RealType logpsi(Psi.evaluateLog(W));
      RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) 
      {
        thisWalker.Age++;
	++nReject; 
        H.rejectedMove(W,thisWalker);
      } 
      else 
      {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
        H.auxHevaluate(W,thisWalker);
	H.saveProperty(thisWalker.getPropertyBase());
	++nAccept;
      }
    }
  }
  
  /// Constructor.
  VMCUpdateAllWithDrift::VMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg) 
    { 
    }

  VMCUpdateAllWithDrift::~VMCUpdateAllWithDrift()
  {
  }

  void VMCUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure) 
  {
    for(;it != it_end;++it) 
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen);
      if(!W.makeMoveWithDrift(thisWalker,deltaR, m_sqrttau)) continue;
      //W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
      //W.update();
      RealType logpsi(Psi.evaluateLog(W));
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      setScaledDrift(m_tauovermass,W.G,drift);

      //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      
      RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) 
      {
        thisWalker.Age++;
	++nReject; 
	H.rejectedMove(W,thisWalker);
      } 
      else 
      {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
	thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
        H.auxHevaluate(W,thisWalker);
	H.saveProperty(thisWalker.getPropertyBase());
	++nAccept;
      }
    }
  }
    /*
  * WFMCUpdateAllWithKill member functions
  */
  /// Constructor.
  WFMCUpdateAllWithReweight::WFMCUpdateAllWithReweight(MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, int length, int index ): QMCUpdateBase(w,psi,h,rg)
    {
      Elength=length-1;
      Eindex =index; 
    }
  
  /// destructor
  WFMCUpdateAllWithReweight::~WFMCUpdateAllWithReweight() { }

  //void DMCUpdateAllWithKill::initWalkers(WalkerIter_t it, WalkerIter_t it_end){
  //}

  //void DMCUpdateAllWithKill::updateWalkers(WalkerIter_t it, WalkerIter_t it_end){
  //}
  /** advance all the walkers with killnode==yes
   */
  void WFMCUpdateAllWithReweight::advanceWalkers(WalkerIter_t it
      , WalkerIter_t it_end, bool measure) 
  {
    for(;it != it_end;++it) 
    {

      Walker_t& thisWalker(**it);
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,RandomGen);
      //reject illegal positions or big displacement
      if(!W.makeMoveWithDrift(thisWalker,deltaR, m_sqrttau)) continue;
      ///W.R = m_sqrttau*deltaR + thisWalker.Drift;
      ///W.R += thisWalker.R;
      ///W.update();
      //save old local energy
      RealType eold = thisWalker.Properties(LOCALENERGY);
      RealType enew = eold;
      //evaluate wave function
      RealType logpsi(Psi.evaluateLog(W));

      bool accepted=false;
      enew=H.evaluate(W);
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);

      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob= std::min(std::exp(logGb-logGf +2.0*(logpsi-thisWalker.Properties(LOGPSI))),1.0);
      if(RandomGen() > prob)
      {
        enew=eold;
        thisWalker.Age++;
        //           thisWalker.Properties(R2ACCEPTED)=0.0;
        //           thisWalker.Properties(R2PROPOSED)=rr_proposed;
        H.rejectedMove(W,thisWalker);
      } 
      else 
      {
        thisWalker.Age=0;
        accepted=true;  
        thisWalker.R = W.R;
        thisWalker.Drift = drift;
        //           rr_accepted = rr_proposed;
        //           thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
        thisWalker.resetProperty(logpsi,Psi.getPhase(),enew);
        H.auxHevaluate(W,thisWalker);
        H.saveProperty(thisWalker.getPropertyBase());
      }
      thisWalker.Weight *= std::exp(-Tau*( enew -thisWalker.PropertyHistory[Eindex][Elength]));
      thisWalker.addPropertyHistoryPoint(Eindex,enew);

      if(accepted) 
        ++nAccept;
      else 
        ++nReject;
    }
  }
  
  
}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
