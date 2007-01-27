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

  VMCUpdateAll::VMCUpdateAll(ParticleSet& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg), nSubSteps(1)
    { 
    }

  VMCUpdateAll::~VMCUpdateAll()
  {
  }

  bool VMCUpdateAll::put(xmlNodePtr cur)
  {
    ParameterSet params;
    params.add(nSubSteps,"subSteps","int"); params.add(nSubSteps,"substeps","int");
    return params.put(cur);
  }

  void VMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {
    while(it != it_end) 
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen);
      W.R = m_sqrttau*deltaR + thisWalker.R;
      W.update();
      RealType logpsi(Psi.evaluateLog(W));
      RealType g= exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) {
        thisWalker.Age++;
	++nReject; 
      } else {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
	H.saveProperty(thisWalker.getPropertyBase());
	++nAccept;
      }
      ++it; 
    }
  }
  
  /// Constructor.
  VMCUpdateAllWithDrift::VMCUpdateAllWithDrift(ParticleSet& w, TrialWaveFunction& psi, 
      QMCHamiltonian& h, RandomGenerator_t& rg):
    QMCUpdateBase(w,psi,h,rg) 
    { 
    }

  VMCUpdateAllWithDrift::~VMCUpdateAllWithDrift()
  {
  }

  void VMCUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end) 
  {
    while(it != it_end) 
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);

      makeGaussRandomWithEngine(deltaR,RandomGen);
      
      W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
      
      W.update();
      
      RealType logpsi(Psi.evaluateLog(W));
 
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      
      setScaledDrift(Tau,W.G,drift);
      
      //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      
      RealType g= exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) {
        thisWalker.Age++;
	++nReject; 
      } else {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
	thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
	H.saveProperty(thisWalker.getPropertyBase());
	++nAccept;
      }
      ++it; 
    }
  }
}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
