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

namespace qmcplusplus
  {

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
    for (; it!= it_end; ++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        bool Cont(true);
        if (!W.makeMove(thisWalker,deltaR, m_sqrttau)) Cont=false;
        if (!Cont)
          {
            H.rejectedMove(W,thisWalker);
            continue;
          }

        //W.R = m_sqrttau*deltaR + thisWalker.R;
        //W.update();

        RealType logpsi(Psi.evaluateLog(W));
//         W.Properties(LOGPSI)=logpsi;
        RealType g= std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
        if (RandomGen() > g)
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
    for (;it != it_end;++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        bool Cont(true);
        if (!W.makeMoveWithDrift(thisWalker,deltaR, m_sqrttau)) Cont=false;
        if (!Cont)
          {
            H.rejectedMove(W,thisWalker);
            continue;
          }
        //W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
        //W.update();
        RealType logpsi(Psi.evaluateLog(W));
        RealType logGf = -0.5*Dot(deltaR,deltaR);
        setScaledDrift(m_tauovermass,W.G,drift);

        //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
        deltaR = thisWalker.R - W.R - drift;
        RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);

        RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
        if (RandomGen() > g)
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



  VMCUpdateAllSamplePsi::VMCUpdateAllSamplePsi(MCWalkerConfiguration& w, TrialWaveFunction& psi,
      QMCHamiltonian& h, RandomGenerator_t& rg):
      QMCUpdateBase(w,psi,h,rg), nSubSteps(1)
  {
    myParams.add(nSubSteps,"subSteps","int");
    myParams.add(nSubSteps,"substeps","int");
  }

  VMCUpdateAllSamplePsi::~VMCUpdateAllSamplePsi()
  {
  }

  void VMCUpdateAllSamplePsi::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
  {
    for (; it!= it_end; ++it)
      {
        MCWalkerConfiguration::Walker_t& thisWalker(**it);
        makeGaussRandomWithEngine(deltaR,RandomGen);
        if (!W.makeMove(thisWalker,deltaR,m_sqrttau)) continue;
        //W.R = m_sqrttau*deltaR + thisWalker.R;
        //W.update();

        RealType logpsi(Psi.evaluateLog(W));
        RealType g= std::exp((logpsi-thisWalker.Properties(LOGPSI)));
        if (RandomGen() > g)
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
            thisWalker.Weight=std::exp(logpsi);
            H.saveProperty(thisWalker.getPropertyBase());
            ++nAccept;
          }
      }
  }

}

/***************************************************************************
 * $RCSfile: VMCUpdateAll.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdateAll.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
