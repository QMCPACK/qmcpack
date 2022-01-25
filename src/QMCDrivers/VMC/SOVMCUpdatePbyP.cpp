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


#include "SOVMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Concurrency/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif


namespace qmcplusplus
{
/// Constructor
SOVMCUpdatePbyP::SOVMCUpdatePbyP(MCWalkerConfiguration& w,
                                 TrialWaveFunction& psi,
                                 QMCHamiltonian& h,
                                 RandomGenerator& rg)
    : QMCUpdateBase(w, psi, h, rg),
      buffer_timer_(*timer_manager.createTimer("SOVMCUpdatePbyP::Buffer", timer_level_medium)),
      movepbyp_timer_(*timer_manager.createTimer("SOVMCUpdatePbyP::MovePbyP", timer_level_medium)),
      hamiltonian_timer_(*timer_manager.createTimer("SOVMCUpdatePbyP::Hamiltonian", timer_level_medium)),
      collectables_timer_(*timer_manager.createTimer("SOVMCUpdatePbyP::Collectables", timer_level_medium))
{}

SOVMCUpdatePbyP::~SOVMCUpdatePbyP() {}

void SOVMCUpdatePbyP::advanceWalker(Walker_t& thisWalker, bool recompute)
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
    makeGaussRandomWithEngine(deltaS, RandomGen);
    moved = false;
    for (int ig = 0; ig < W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau * MassInvS[ig];
      RealType oneover2tau = 0.5 / (tauovermass);
      RealType sqrttau     = std::sqrt(tauovermass);
      Psi.prepareGroup(W, ig);
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        PosType dr;
        ParticleSet::Scalar_t ds;
        if (UseDrift)
        {
          ComplexType spingrad_now;
          GradType grad_now = Psi.evalGradWithSpin(W, iat, spingrad_now);
          DriftModifier->getDrift(tauovermass, grad_now, dr);
          DriftModifier->getDrift(tauovermass / spinMass, spingrad_now, ds);
          dr += sqrttau * deltaR[iat];
          ds += std::sqrt(tauovermass / spinMass) * deltaS[iat];
        }
        else
        {
          dr = sqrttau * deltaR[iat];
          ds = std::sqrt(tauovermass / spinMass) * deltaS[iat];
        }
        if (!W.makeMoveAndCheckWithSpin(iat, dr, ds))
        {
          ++nReject;
          W.accept_rejectMove(iat, false);
          continue;
        }
        RealType prob(0);
        if (UseDrift)
        {
          GradType grad_new;
          ComplexType spingrad_new;
          prob = std::norm(Psi.calcRatioGradWithSpin(W, iat, grad_new, spingrad_new));
          DriftModifier->getDrift(tauovermass, grad_new, dr);
          dr             = W.R[iat] - W.getActivePos() - dr;
          RealType logGb = -oneover2tau * dot(dr, dr);
          RealType logGf = mhalf * dot(deltaR[iat], deltaR[iat]);

          DriftModifier->getDrift(tauovermass / spinMass, spingrad_new, ds);
          ds = W.spins[iat] - W.getActiveSpinVal() - ds;
          logGb += -spinMass * oneover2tau * ds * ds;
          logGf += mhalf * deltaS[iat] * deltaS[iat];

          prob *= std::exp(logGb - logGf);
        }
        else
        {
          prob = std::norm(Psi.calcRatio(W, iat));
        }

        bool is_accepted = false;
        if (prob >= std::numeric_limits<RealType>::epsilon() && RandomGen() < prob)
        {
          is_accepted = true;
          moved       = true;
          ++nAccept;
          Psi.acceptMove(W, iat, true);
        }
        else
        {
          ++nReject;
          Psi.rejectMove(iat);
        }
        W.accept_rejectMove(iat, is_accepted);
      }
    }
    Psi.completeUpdates();
  }
  W.donePbyP();
  movepbyp_timer_.stop();
  buffer_timer_.start();
  RealType logpsi = Psi.updateBuffer(W, w_buffer, recompute);
  if (debug_checks_ & DriverDebugChecks::CHECKGL_AFTER_MOVES)
    checkLogAndGL(W, Psi, "checkGL_after_moves");
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
