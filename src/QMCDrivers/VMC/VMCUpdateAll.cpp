//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

VMCUpdateAll::VMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                           QMCHamiltonian& h, RandomGenerator_t& rg)
  : QMCUpdateBase(w,psi,h,rg)
{
  UpdatePbyP=false;
}

VMCUpdateAll::~VMCUpdateAll()
{
}

void VMCUpdateAll::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  /* thisWalker.R will track the last accepted configuration
   * W.R will track the proposed configuration 
   *
   * upon Call:
   *   thisWalker.R,G,L,DistTables,SK,Properties must be consistent
   *   recompute flag is for pbyp driver and is not used here
   *
   * upon Return:
   *   thisWalker.R,G,L,DistTables,SK,Properties must be kept consistent
   * */
  bool updated=false;
  W.loadWalker(thisWalker,false); // W.R,G,L = thisWalker.R,G,L
  RealType logpsi_old=thisWalker.Properties(LOGPSI);
  for (int iter=0; iter<nSubSteps; ++iter)
  { // make a few Monte-Carlo steps to decorrelate samples without calculating observables
    makeGaussRandomWithEngine(deltaR,RandomGen); // fill deltaR
    if (W.makeMove(thisWalker,deltaR,SqrtTauOverMass))
    { // W.R += dR*dt; W.DistTables,SK are updated; W.G,L are now stale
      RealType logpsi=Psi.evaluateLog(W); // update W.G,L at changed W.R; update Psi.LogValue,PhaseValue
      RealType g= std::exp(2.0*(logpsi-logpsi_old));
      if (RandomGen() > g)
      { // logpsi_old is still good
        thisWalker.Age++; ++nReject; updated=false;
      }
      else
      { // move is accepted; logpsi_old and R_old are stale
        logpsi_old=logpsi; // update logpsi_old
        thisWalker.R=W.R;  // update R_old; side effect: thisWalker.G,L,DistTables,SK,Properties are now stale
        ++nAccept; updated=true;
      }
    }
    else // if a proposed move is too large, it will be pre-emptively rejected
    { // W.R,G,L,DistTables,SK are not updated. i.e. they are still consistent
      thisWalker.Age++; ++nReject; updated=false;
    }
  }

  if(!updated)
  { // W.G and W.L have to be computed because the last move was rejected
    W.update(thisWalker.R); // move W back to last accepted configuration; W.DistTables, SK are updated
    logpsi_old=Psi.evaluateLog(W); // update W.G,L
  } // W and logpsi_old are up-to-date at this point

  RealType eloc = H.evaluate(W); // calculate local energy; W.SK must be up-to-date if Coulomb interaction is used with periodic boundary. W.SK is used to calculate the long-range part of the Coulomb potential.
  thisWalker.R = W.R;
  thisWalker.resetProperty(logpsi_old,Psi.getPhase(),eloc); // update thisWalker::Properties[LOGPSI,SIGN,LOCALENERGY]
  H.auxHevaluate(W,thisWalker); // update auxiliary observables, i.e. fill H::Observables
  H.saveProperty(thisWalker.getPropertyBase()); // copy H::Observables to thisWalker::Properties
}

/// Constructor.
VMCUpdateAllWithDrift::VMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg)
{
  UpdatePbyP=false;
}

VMCUpdateAllWithDrift::~VMCUpdateAllWithDrift()
{
}

void VMCUpdateAllWithDrift::advanceWalker(Walker_t& thisWalker, bool recompute)
{ /* see VMCUpdateAll::VMCUpdateAll */

  bool updated=false;
  W.loadWalker(thisWalker,false);
  RealType logpsi_old=thisWalker.Properties(LOGPSI);
  // R_old will be stored in thisWalker.R

  for (int iter=0; iter<nSubSteps; ++iter)
  { // make a few Monte-Carlo steps to decorrelate samples without calculating observables
    assignDrift(Tau,MassInvP,W.G,drift); // fill variable drift, require W.G to be up-to-date
    makeGaussRandomWithEngine(deltaR,RandomGen); // fill variable deltaR
    if (W.makeMoveWithDrift(thisWalker,drift,deltaR,SqrtTauOverMass))
    { // W.R = thisWalker.R + drift + deltaR; W.DistTables,SK are updated; W.G,L are now stale
      
      RealType logpsi=Psi.evaluateLog(W);  // update W.G,L; update Psi.PhaseValue,LogValue
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      assignDrift(Tau,MassInvP,W.G,drift); // update drift at proposed configuration
      deltaR = thisWalker.R - W.R - drift; // hijack deltaR to hold reverse move
      RealType logGb=logBackwardGF(deltaR);

      RealType g= std::exp(logGb-logGf+2.0*(logpsi-logpsi_old));
      // accept or reject
      if (RandomGen() > g)
      {
        W.G=thisWalker.G; // update W.G to last accepted configuration
        thisWalker.Age++; ++nReject; updated=false;
      }
      else
      {
        thisWalker.R = W.R; thisWalker.G=W.G;
        ++nAccept; logpsi_old=logpsi; updated=true;
      }
    }
    else
    { // if a proposed move is too large, it will be pre-emptively rejected
      W.G=thisWalker.G; // is this necessary?
      thisWalker.Age++; ++nReject; updated=false;
    }
  }

  if(!updated)
  { // W.G and W.L have to be computed because the last move was rejected
    W.update(thisWalker.R);
    logpsi_old=Psi.evaluateLog(W);
  }

  RealType eloc = H.evaluate(W); // calculate local energy; W.SK must be up-to-date if Coulomb interaction is used with periodic boundary. W.SK is used to calculate the long-range part of the Coulomb potential.
  thisWalker.R = W.R;
  thisWalker.resetProperty(logpsi_old,Psi.getPhase(),eloc);
  H.auxHevaluate(W,thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
}

VMCUpdateAllWithDrift::RealType VMCUpdateAllWithDrift::advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios)
{
  // !!!! 2017-10-18: fixed SK update bug in advanceWalker. This function needs to be checked!
  W.loadWalker(w1,false);
  std::vector<RealType> orb_ratios(4,1.0);
  int nshells(rs.size());
  std::vector<RealType> orb_phases0,orb_logs0;
  Psi.evaluateLog(W);
  RealType pt=Psi.getPhase();
  RealType lt=Psi.getLogPsi();
  Psi.getPhases(orb_phases0);
  Psi.getLogs(orb_logs0);
  int nmoved(0);
  int i(0);
  while(i<rs.size())
  {
    deltaR=0;
    while(nmoved<rs[i])
    {
      int iat=iats[nmoved];
      W.R[iat]+=dR[nmoved];
      nmoved++;
    }
    RealType logpsi(Psi.evaluateLog(W));
    std::vector<RealType> orb_phases,orb_logs;
    Psi.getPhases(orb_phases);
    Psi.getLogs(orb_logs);
    orb_ratios[0]=std::cos(pt-Psi.getPhase())*std::exp(logpsi-lt);
    orb_ratios[1]=std::cos(orb_phases[0] - orb_phases0[0])*std::exp(orb_logs[0] - orb_logs0[0]);
    orb_ratios[2]=std::cos(orb_phases[1] - orb_phases0[1])*std::exp(orb_logs[1] - orb_logs0[1]);
    orb_ratios[3]=orb_ratios[0]/orb_ratios[1];
    while(nmoved==rs[i])
    {
      ratios[i]=orb_ratios[0];
      ratios[nshells+i]=orb_ratios[1];
      ratios[nshells*2+i]=orb_ratios[2];
      ratios[nshells*3+i]=orb_ratios[3];
      i++;
      if(i==rs.size())
        break;
    }
  }
  return 0.0;
}

}
