//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCUpdateBase.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

/** advance a walker: walker move, use drift and vmc
 */
void QMCUpdateBase::advanceWalker(Walker_t& thisWalker)
{
  W.loadWalker(thisWalker,false);
  RealType logpsi0=thisWalker.Properties(LOGPSI);
  RealType phase0=thisWalker.Properties(SIGN);
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    if (!W.makeMoveWithDrift(thisWalker, drift, deltaR, m_sqrttau))
    {
      ++nReject;
      continue;
    }
    RealType logpsi(Psi.evaluateLog(W));
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    // setScaledDrift(m_tauovermass,W.G,drift);
    nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
    deltaR = thisWalker.R - W.R - drift;
    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
    RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
    }
    else
    {
      thisWalker.R=W.R;
      thisWalker.G=W.G;
      thisWalker.L=W.L;
      //skip energy
      thisWalker.resetProperty(logpsi,Psi.getPhase(),0);
      //update logpsi0, phase0
      logpsi0=logpsi;
      phase0=Psi.getPhase();
      ++nAccept;
    }
  }
  //measure energy
  W.loadWalker(thisWalker,true);
  RealType eloc=H.evaluate(W);
  thisWalker.resetProperty(logpsi0,phase0,eloc);
  H.auxHevaluate(W,thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
}

/** advance of a walker using VMC+drift */
void QMCUpdateBase::advancePbyP(Walker_t& thisWalker)
{
  Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
  W.loadWalker(thisWalker,true);
  Psi.copyFromBuffer(W,w_buffer);
  bool moved = false;
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    for(int ig=0; ig<W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau*MassInvS[ig];
      RealType oneover2tau = 0.5/(tauovermass);
      RealType sqrttau = std::sqrt(tauovermass);
      for (int iat=W.first(ig); iat<W.last(ig); ++iat)
      {
        GradType grad_now=Psi.evalGrad(W,iat), grad_new;
        mPosType dr;
        getScaledDrift(tauovermass,grad_now,dr);
        dr += sqrttau*deltaR[iat];
        if (!W.makeMoveAndCheck(iat,dr))
        {
          ++nReject;
          continue;
        }
        //PosType newpos = W.makeMove(iat,dr);
        RealType ratio = Psi.ratioGrad(W,iat,grad_new);
        RealType prob = ratio*ratio;
        //zero is always rejected
        if (prob<std::numeric_limits<RealType>::epsilon())
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
          continue;
        }
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
        getScaledDrift(tauovermass,grad_new,dr);
        dr = thisWalker.R[iat]-W.R[iat]-dr;
        RealType logGb = -oneover2tau*dot(dr,dr);
        //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
        if (RandomGen() < prob*std::exp(logGb-logGf))
        {
          moved = true;
          ++nAccept;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
        }
        else
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
    //for subSteps must update thiswalker
    thisWalker.R=W.R;
    thisWalker.G=W.G;
    thisWalker.L=W.L;
  }
  //Always compute the energy
  {
    RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
    W.saveWalker(thisWalker);
    RealType eloc=H.evaluate(W);
    //thisWalker.resetProperty(std::log(std::abs(psi)), psi,eloc);
    thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
    H.auxHevaluate(W,thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());
  }
  if(!moved)
    ++nAllRejected;
}
}

