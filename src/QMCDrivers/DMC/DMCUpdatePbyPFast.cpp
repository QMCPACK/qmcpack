//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif
//#define TEST_INNERBRANCH


namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

TimerNameList_t<DMCTimers> DMCTimerNames = {{DMC_buffer, "DMCUpdatePbyP::Buffer"},
                                            {DMC_movePbyP, "DMCUpdatePbyP::movePbyP"},
                                            {DMC_hamiltonian, "DMCUpdatePbyP::Hamiltonian"},
                                            {DMC_collectables, "DMCUpdatePbyP::Collectables"},
                                            {DMC_tmoves, "DMCUpdatePbyP::Tmoves"}};


/// Constructor.
DMCUpdatePbyPWithRejectionFast::DMCUpdatePbyPWithRejectionFast(MCWalkerConfiguration& w,
                                                               TrialWaveFunction& psi,
                                                               QMCHamiltonian& h,
                                                               RandomGenerator_t& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  setup_timers(myTimers, DMCTimerNames, timer_level_medium);
}

/// destructor
DMCUpdatePbyPWithRejectionFast::~DMCUpdatePbyPWithRejectionFast() {}

void DMCUpdatePbyPWithRejectionFast::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  {
    ScopedTimer local_timer(myTimers[DMC_buffer]);
    W.loadWalker(thisWalker, true);
    Psi.copyFromBuffer(W, w_buffer);
  }
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  FullPrecRealType eold(thisWalker.Properties(WP::LOCALENERGY));
  FullPrecRealType enew(eold);
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  RealType gf_acc      = 1.0;
  {
    ScopedTimer local_timer(myTimers[DMC_movePbyP]);
    for (int ig = 0; ig < W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau * MassInvS[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      Psi.prepareGroup(W, ig);
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        //get the displacement
        GradType grad_iat = Psi.evalGrad(W, iat);
        PosType dr;
        DriftModifier->getDrift(tauovermass, grad_iat, dr);
        dr += sqrttau * deltaR[iat];
        bool is_valid = W.makeMoveAndCheck(iat, dr);
        RealType rr   = tauovermass * dot(deltaR[iat], deltaR[iat]);
        rr_proposed += rr;
        if (!is_valid || rr > m_r2max)
        {
          ++nRejectTemp;
          W.accept_rejectMove(iat, false);
          continue;
        }
        ValueType ratio = Psi.calcRatioGrad(W, iat, grad_iat);
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
          //Use the force of the particle iat
          DriftModifier->getDrift(tauovermass, grad_iat, dr);
          dr                     = W.R[iat] - W.activePos - dr;
          FullPrecRealType logGb = -oneover2tau * dot(dr, dr);
          RealType prob          = std::norm(ratio) * std::exp(logGb - logGf);
          bool is_accepted       = false;
          if (RandomGen() < prob)
          {
            is_accepted = true;
            ++nAcceptTemp;
            Psi.acceptMove(W, iat, true);
            rr_accepted += rr;
            gf_acc *= prob; //accumulate the ratio
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
      ScopedTimer local_timer(myTimers[DMC_buffer]);
      thisWalker.Age = 0;
      logpsi         = Psi.updateBuffer(W, w_buffer, recompute);
      assert(checkLogAndGL(W, Psi));
      W.saveWalker(thisWalker);
    }
    {
      ScopedTimer local_timer(myTimers[DMC_hamiltonian]);
      enew = H.evaluateWithToperator(W);
    }
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 1.0);
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    {
      ScopedTimer local_timer(myTimers[DMC_collectables]);
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
    enew   = eold; //copy back old energy
    gf_acc = 1.0;
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  }
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  {
    ScopedTimer local_timer(myTimers[DMC_tmoves]);
    const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
    if (NonLocalMoveAcceptedTemp > 0)
    {
      RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
      // debugging lines
      //W.update(true);
      //RealType logpsi2 = Psi.evaluateLog(W);
      //if(logpsi!=logpsi2) std::cout << " logpsi " << logpsi << " logps2i " << logpsi2 << " diff " << logpsi2-logpsi << std::endl;
      W.saveWalker(thisWalker);
      NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
    }
  }
  nAccept += nAcceptTemp;
  nReject += nRejectTemp;

  setMultiplicity(thisWalker);
}

} // namespace qmcplusplus
