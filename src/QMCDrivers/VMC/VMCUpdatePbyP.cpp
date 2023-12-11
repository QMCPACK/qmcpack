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


#include "VMCUpdatePbyP.h"
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
VMCUpdatePbyP::VMCUpdatePbyP(MCWalkerConfiguration& w,
                             TrialWaveFunction& psi,
                             QMCHamiltonian& h,
                             RandomBase<FullPrecRealType>& rg)
    : QMCUpdateBase(w, psi, h, rg),
      buffer_timer_(createGlobalTimer("VMCUpdatePbyP::Buffer", timer_level_medium)),
      movepbyp_timer_(createGlobalTimer("VMCUpdatePbyP::MovePbyP", timer_level_medium)),
      hamiltonian_timer_(createGlobalTimer("VMCUpdatePbyP::Hamiltonian", timer_level_medium)),
      collectables_timer_(createGlobalTimer("VMCUpdatePbyP::Collectables", timer_level_medium))
{}

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
      Psi.prepareGroup(W, ig);
      for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      {
        PosType dr;
        if (UseDrift)
        {
          GradType grad_now = Psi.evalGrad(W, iat);
          DriftModifier->getDrift(tauovermass, grad_now, dr);
          dr += sqrttau * deltaR[iat];
        }
        else
          dr = sqrttau * deltaR[iat];

        if (!W.makeMoveAndCheck(iat, dr))
        {
          ++nReject;
          W.accept_rejectMove(iat, false);
          continue;
        }

        RealType prob(0);
        if (UseDrift)
        {
          GradType grad_new;
          prob = std::norm(Psi.calcRatioGrad(W, iat, grad_new));
          DriftModifier->getDrift(tauovermass, grad_new, dr);
          dr             = W.R[iat] - W.getActivePos() - dr;
          RealType logGb = -oneover2tau * dot(dr, dr);
          RealType logGf = mhalf * dot(deltaR[iat], deltaR[iat]);
          prob *= std::exp(logGb - logGf);
        }
        else
          prob = std::norm(Psi.calcRatio(W, iat));

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
