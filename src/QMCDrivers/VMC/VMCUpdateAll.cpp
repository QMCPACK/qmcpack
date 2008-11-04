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
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    measure &= (compEstimator != 0);
    if(measure) compEstimator->startAccumulate();
#endif
    while(it != it_end) 
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);
      makeGaussRandomWithEngine(deltaR,RandomGen);
      W.R = m_sqrttau*deltaR + thisWalker.R;
      W.update();
      RealType logpsi(Psi.evaluateLog(W));
      RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) {
        thisWalker.Age++;
	++nReject; 
        thisWalker.rejectedMove();
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
        if(measure)
        {//evaluate the old value
          W.R = thisWalker.R;
          W.update();
          compEstimator->accumulate(W,1.0);
          
        }
#endif
      } else {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
        H.auxHevaluate(W,thisWalker);
	H.saveProperty(thisWalker.getPropertyBase());
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
        if(measure) compEstimator->accumulate(W,1.0);
#endif
	++nAccept;
      }
      ++it; 
    }
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(measure) compEstimator->stopAccumulate(-1);
#endif
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
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    measure &= (compEstimator != 0);
    if(measure) compEstimator->startAccumulate();
#endif

    while(it != it_end) 
    {
      MCWalkerConfiguration::Walker_t& thisWalker(**it);

      makeGaussRandomWithEngine(deltaR,RandomGen);
      
      W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
      
      W.update();
      
      RealType logpsi(Psi.evaluateLog(W));
 
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      
      setScaledDrift(m_tauovermass,W.G,drift);
      
      //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      
      RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
      if(RandomGen() > g) {
        thisWalker.Age++;
	++nReject; 
        thisWalker.rejectedMove();
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
        if(measure)
        {//evaluate the old value
          W.R = thisWalker.R;
          W.update();
          compEstimator->accumulate(W,1.0);
        }
#endif
      } else {
        RealType eloc=H.evaluate(W);
	thisWalker.R = W.R;
	thisWalker.Drift = drift;
        thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
        H.auxHevaluate(W,thisWalker);
	H.saveProperty(thisWalker.getPropertyBase());
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
        if(measure) compEstimator->accumulate(W,1.0);
#endif
	++nAccept;
      }
      ++it; 
    }
#if defined(ENABLE_COMPOSITE_ESTIMATOR)
    if(measure) compEstimator->stopAccumulate(-1);
#endif
  }
}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $ 
 ***************************************************************************/
