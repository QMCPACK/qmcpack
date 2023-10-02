//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "DMCUpdateAll.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/WalkerProperties.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/// Constructor.
DMCUpdateAllWithRejection::DMCUpdateAllWithRejection(MCWalkerConfiguration& w,
                                                     TrialWaveFunction& psi,
                                                     QMCHamiltonian& h,
                                                     RandomBase<FullPrecRealType>& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  UpdatePbyP = false;
}

/// destructor
DMCUpdateAllWithRejection::~DMCUpdateAllWithRejection() {}

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void DMCUpdateAllWithRejection::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  W.loadWalker(thisWalker, false);
  //create a 3N-Dimensional Gaussian with variance=1
  RealType nodecorr = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
  //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  makeGaussRandomWithEngine(deltaR, RandomGen);
  //save old local energy
  RealType eold        = thisWalker.Properties(WP::LOCALENERGY);
  RealType enew        = eold;
  bool accepted        = false;
  RealType rr_accepted = 0.0;
  RealType rr_proposed = 0.0;
  RealType logpsi;

  if (W.makeMoveAllParticlesWithDrift(thisWalker, drift, deltaR, SqrtTauOverMass))
  {
    //evaluate the new wave function
    logpsi = Psi.evaluateLog(W);
    //fixed node
    if (!branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      RealType logGf = -0.5 * Dot(deltaR, deltaR);
      nodecorr       = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
      deltaR         = thisWalker.R - W.R - drift;
      RealType logGb = logBackwardGF(deltaR);
      //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob = std::min(std::exp(logGb - logGf + 2.0 * (logpsi - thisWalker.Properties(WP::LOGPSI))), 1.0);
      //calculate rr_proposed here
      deltaR      = W.R - thisWalker.R;
      rr_proposed = Dot(deltaR, deltaR);
      if (RandomGen() <= prob)
      {
        accepted    = true;
        rr_accepted = rr_proposed;
      }
    }
  }

  // recompute Psi if the move is rejected
  if (!accepted)
  {
    W.loadWalker(thisWalker, false);
    W.update();
    logpsi = Psi.evaluateLog(W);
  }

  // evaluate Hamiltonian
  enew = H.evaluateWithToperator(W);
  H.auxHevaluate(W, thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());

  // operate on thisWalker.
  if (accepted)
  {
    W.saveWalker(thisWalker);
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, nodecorr);
  }
  else
  {
    thisWalker.Age++;
    thisWalker.Properties(WP::R2ACCEPTED) = 0.0;
    thisWalker.Properties(WP::R2PROPOSED) = rr_proposed;
  }

  const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
  if (NonLocalMoveAcceptedTemp > 0)
  {
    W.saveWalker(thisWalker);
    thisWalker.resetProperty(Psi.getLogPsi(), Psi.getPhase(), enew);
    // debugging lines
    //logpsi = Psi.getLogPsi();
    //W.update(true);
    //RealType logpsi2 = Psi.evaluateLog(W);
    //if(logpsi!=logpsi2) std::cout << " logpsi " << logpsi << " logpsi2 " << logpsi2
    //                              << " diff " << logpsi2-logpsi << std::endl;

    NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
  }

  thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  //branchEngine->accumulate(eold,1);
  if (accepted)
    ++nAccept;
  else
    ++nReject;

  setMultiplicity(thisWalker);
}

/*
 * DMCUpdateAllWithKill member functions
 */
/// Constructor.
DMCUpdateAllWithKill::DMCUpdateAllWithKill(MCWalkerConfiguration& w,
                                           TrialWaveFunction& psi,
                                           QMCHamiltonian& h,
                                           RandomBase<FullPrecRealType>& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  UpdatePbyP = false;
}

/// destructor
DMCUpdateAllWithKill::~DMCUpdateAllWithKill() {}

/** advance all the walkers with killnode==yes
 */
void DMCUpdateAllWithKill::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  W.loadWalker(thisWalker, false);
  //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  RealType nodecorr = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  //if(!W.makeMoveAllParticlesWithDrift(thisWalker,drift,deltaR, m_sqrttau))
  if (!W.makeMoveAllParticlesWithDrift(thisWalker, drift, deltaR, SqrtTauOverMass))
  {
    H.rejectedMove(W, thisWalker);
    return;
  }
  //save old local energy
  RealType eold = thisWalker.Properties(WP::LOCALENERGY);
  RealType enew = eold;
  RealType logpsi(Psi.evaluateLog(W));
  bool accepted        = false;
  RealType rr_accepted = 0.0;
  nodecorr             = 0.0;
  if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
  {
    thisWalker.Age++;
    thisWalker.willDie();
  }
  else
  {
    enew           = H.evaluate(W);
    RealType logGf = -0.5 * Dot(deltaR, deltaR);
    //nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
    nodecorr       = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
    deltaR         = thisWalker.R - W.R - drift;
    RealType logGb = logBackwardGF(deltaR);
    //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
    RealType prob = std::min(std::exp(logGb - logGf + 2.0 * (logpsi - thisWalker.Properties(WP::LOGPSI))), 1.0);
    //calculate rr_proposed here
    deltaR               = W.R - thisWalker.R;
    RealType rr_proposed = Dot(deltaR, deltaR);
    if (RandomGen() > prob)
    {
      enew = eold;
      thisWalker.Age++;
      thisWalker.Properties(WP::R2ACCEPTED) = 0.0;
      thisWalker.Properties(WP::R2PROPOSED) = rr_proposed;
      H.rejectedMove(W, thisWalker);
    }
    else
    {
      thisWalker.Age = 0;
      accepted       = true;
      W.saveWalker(thisWalker);
      rr_accepted = rr_proposed;
      thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, nodecorr);
      H.auxHevaluate(W, thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    //         std::cout <<logpsi<<"  "<<Psi.getPhase()<<"  "<<enew<<"  "<<rr_accepted<<"  "<<rr_proposed<<"  "<<nodecorr<< std::endl;
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  }
  if (accepted)
    ++nAccept;
  else
    ++nReject;

  setMultiplicity(thisWalker);
}
} // namespace qmcplusplus
