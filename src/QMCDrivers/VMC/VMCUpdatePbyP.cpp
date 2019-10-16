//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif


namespace qmcplusplus
{

/// Constructor
VMCUpdatePbyP::VMCUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
    : QMCUpdateBase(w, psi, h, rg),
      buffer_timer_(*TimerManager.createTimer("VMCUpdatePbyP::Buffer",timer_level_medium)),
      movepbyp_timer_(*TimerManager.createTimer("VMCUpdatePbyP::MovePbyP",timer_level_medium)),
      hamiltonian_timer_(*TimerManager.createTimer("VMCUpdatePbyP::Hamiltonian",timer_level_medium)),
      collectables_timer_(*TimerManager.createTimer("VMCUpdatePbyP::Collectables",timer_level_medium))
{
}

VMCUpdatePbyP::~VMCUpdatePbyP() {}

void VMCUpdatePbyP::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  buffer_timer_.start();
  W.loadWalker(thisWalker, true);
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  Psi.copyFromBuffer(W, w_buffer);
  buffer_timer_.stop();

  // start PbyP moves
  movepbyp_timer_.start();
  bool moved = false;
  constexpr RealType mhalf(-0.5);
  for (int iter = 0; iter < nSubSteps; ++iter)
  {
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR, RandomGen);
    moved = false;
    for (int ig = 0; ig < W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau * MassInvS[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        W.setActive(iat);
        PosType dr;
        if (UseDrift)
        {
          GradType grad_now = Psi.evalGrad(W, iat);
          DriftModifier->getDrift(tauovermass, grad_now, dr);
          dr += sqrttau * deltaR[iat];
        }
        else
        {
          dr = sqrttau * deltaR[iat];
        }
        if (!W.makeMoveAndCheck(iat, dr))
        {
          ++nReject;
          continue;
        }
        RealType logGf(1), logGb(1), prob;
        if (UseDrift)
        {
          GradType grad_new;
          RealType ratio = Psi.ratioGrad(W, iat, grad_new);
          prob           = ratio * ratio;
          logGf          = mhalf * dot(deltaR[iat], deltaR[iat]);
          DriftModifier->getDrift(tauovermass, grad_new, dr);
          dr    = W.R[iat] - W.activePos - dr;
          logGb = -oneover2tau * dot(dr, dr);
        }
        else
        {
          RealType ratio = Psi.ratio(W, iat);
          prob           = ratio * ratio;
        }
        if (prob >= std::numeric_limits<RealType>::epsilon() && RandomGen() < prob * std::exp(logGb - logGf))
        {
          moved = true;
          ++nAccept;
          Psi.acceptMove(W, iat);
          W.acceptMove(iat);
        }
        else
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
    Psi.completeUpdates();
  }
  W.donePbyP();
  movepbyp_timer_.stop();
  buffer_timer_.start();
  RealType logpsi = Psi.updateBuffer(W, w_buffer, recompute);
  W.saveWalker(thisWalker);
  buffer_timer_.stop();
  // end PbyP moves
  hamiltonian_timer_.start();
  FullPrecRealType eloc = H.evaluate(W);
  thisWalker.resetProperty(logpsi, Psi.getPhase(), eloc);
  hamiltonian_timer_.stop();
  collectables_timer_.start();
  H.auxHevaluate(W, thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
  collectables_timer_.stop();
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  if (!moved)
    ++nAllRejected;
}

} // namespace qmcplusplus
