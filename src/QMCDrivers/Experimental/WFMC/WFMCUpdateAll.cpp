//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/WFMC/WFMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/// Constructor.
WFMCUpdateAllWithReweight::WFMCUpdateAllWithReweight(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg, int length, int index): QMCUpdateBase(w,psi,h,rg)
{
  Elength=length-1;
  Eindex =index;
}

/// destructor
WFMCUpdateAllWithReweight::~WFMCUpdateAllWithReweight() { }

/** advance all the walkers with killnode==yes
 */
void WFMCUpdateAllWithReweight::advanceWalkers(WalkerIter_t it
    , WalkerIter_t it_end, bool measure)
{
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,false);
    setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    //reject illegal positions or big displacement
    if (!W.makeMoveWithDrift(thisWalker,drift,deltaR, m_sqrttau))
    {
      H.rejectedMove(W,thisWalker);
      H.auxHevaluate(W,thisWalker);
      continue;
    }
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
    if (RandomGen() > prob)
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
      W.saveWalker(thisWalker);
      //           rr_accepted = rr_proposed;
      //           thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    thisWalker.Weight *= std::exp(-Tau*(enew -thisWalker.PropertyHistory[Eindex][Elength]));
    thisWalker.addPropertyHistoryPoint(Eindex,enew);
    if (accepted)
      ++nAccept;
    else
      ++nReject;
  }
}


}

