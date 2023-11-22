//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/DMC/SODMCUpdatePbyP.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCDrivers/DriftOperators.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif
//#define TEST_INNERBRANCH


namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

const TimerNameList_t<SODMCTimers> SODMCTimerNames = {{SODMC_buffer, "SODMCUpdatePbyP::Buffer"},
                                                      {SODMC_movePbyP, "SODMCUpdatePbyP::movePbyP"},
                                                      {SODMC_hamiltonian, "SODMCUpdatePbyP::Hamiltonian"},
                                                      {SODMC_collectables, "SODMCUpdatePbyP::Collectables"},
                                                      {SODMC_tmoves, "SODMCUpdatePbyP::Tmoves"}};


/// Constructor.
SODMCUpdatePbyPWithRejectionFast::SODMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w,
                                                                   TrialWaveFunction& psi,
                                                                   QMCHamiltonian& h,
                                                                   RandomBase<FullPrecRealType>& rg)
    : QMCUpdateBase(w, psi, h, rg), myTimers(getGlobalTimerManager(), SODMCTimerNames, timer_level_medium)

{}

/// destructor
SODMCUpdatePbyPWithRejectionFast::~SODMCUpdatePbyPWithRejectionFast() {}

void SODMCUpdatePbyPWithRejectionFast::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  {
    ScopedTimer local_timer(myTimers[SODMC_buffer]);
    W.loadWalker(thisWalker, true);
    Psi.copyFromBuffer(W, w_buffer);
  }
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  makeGaussRandomWithEngine(deltaS, RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  FullPrecRealType eold(thisWalker.Properties(WP::LOCALENERGY));
  FullPrecRealType enew(eold);
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  {
    ScopedTimer local_timer(myTimers[SODMC_movePbyP]);
    for (int ig = 0; ig < W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau * MassInvS[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      Psi.prepareGroup(W, ig);
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        //get the displacement
        ComplexType spingrad_iat;
        GradType grad_iat = Psi.evalGradWithSpin(W, iat, spingrad_iat);
        PosType dr;
        ParticleSet::Scalar_t ds;
        DriftModifier->getDrift(tauovermass, grad_iat, dr);
        DriftModifier->getDrift(tauovermass / spinMass, spingrad_iat, ds);
        dr += sqrttau * deltaR[iat];
        ds += std::sqrt(tauovermass / spinMass) * deltaS[iat];
        bool is_valid = W.makeMoveAndCheckWithSpin(iat, dr, ds);
        RealType rr   = tauovermass * dot(deltaR[iat], deltaR[iat]);
        rr_proposed += rr;
        if (!is_valid || rr > m_r2max)
        {
          ++nRejectTemp;
          W.accept_rejectMove(iat, false);
          continue;
        }
        ValueType ratio = Psi.calcRatioGradWithSpin(W, iat, grad_iat, spingrad_iat);
        //node is crossed reject the move
        if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.accept_rejectMove(iat, false);
          Psi.rejectMove(iat);
        }
        else
        {
          FullPrecRealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
          logGf += -0.5 * deltaS[iat] * deltaS[iat];
          //Use the force of the particle iat
          DriftModifier->getDrift(tauovermass, grad_iat, dr);
          DriftModifier->getDrift(tauovermass / spinMass, spingrad_iat, ds);
          dr                     = W.R[iat] - W.getActivePos() - dr;
          ds                     = W.spins[iat] - W.getActiveSpinVal() - ds;
          FullPrecRealType logGb = -oneover2tau * dot(dr, dr);
          logGb += -spinMass * oneover2tau * ds * ds;
          RealType prob    = std::norm(ratio) * std::exp(logGb - logGf);
          bool is_accepted = false;

          if (RandomGen() < prob)
          {
            is_accepted = true;

            ++nAcceptTemp;
            Psi.acceptMove(W, iat, true);
            rr_accepted += rr;
          }
          else
          {
            ++nRejectTemp;
            Psi.rejectMove(iat);
          }
          W.accept_rejectMove(iat, is_accepted);
        }
      }
    }
    Psi.completeUpdates();
    W.donePbyP();
  }

  if (nAcceptTemp > 0)
  {
    //need to overwrite the walker properties
    RealType logpsi(0);
    {
      ScopedTimer local_timer(myTimers[SODMC_buffer]);
      thisWalker.Age = 0;
      logpsi         = Psi.updateBuffer(W, w_buffer, recompute);
      if (debug_checks_ & DriverDebugChecks::CHECKGL_AFTER_MOVES)
        checkLogAndGL(W, Psi, "checkGL_after_moves");
      W.saveWalker(thisWalker);
    }
    {
      ScopedTimer local_timer(myTimers[SODMC_hamiltonian]);
      enew = H.evaluateWithToperator(W);
    }
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 1.0);
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    {
      ScopedTimer local_timer(myTimers[SODMC_collectables]);
      H.auxHevaluate(W, thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
  }
  else
  {
    //all moves are rejected: does not happen normally with reasonable wavefunctions
    thisWalker.Age++;
    thisWalker.Properties(WP::R2ACCEPTED) = 0.0;
    //weight is set to 0 for traces
    // consistent w/ no evaluate/auxHevaluate
    RealType wtmp     = thisWalker.Weight;
    thisWalker.Weight = 0.0;
    H.rejectedMove(W, thisWalker);
    thisWalker.Weight = wtmp;
    ++nAllRejected;
    enew = eold; //copy back old energy
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  }
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  {
    ScopedTimer local_timer(myTimers[SODMC_tmoves]);
    const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
    if (NonLocalMoveAcceptedTemp > 0)
    {
      RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
      W.saveWalker(thisWalker);
      NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
    }
  }
  nAccept += nAcceptTemp;
  nReject += nRejectTemp;

  setMultiplicity(thisWalker);
}

} // namespace qmcplusplus
