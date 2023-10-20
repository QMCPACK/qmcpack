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


#include "SOVMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

SOVMCUpdateAll::SOVMCUpdateAll(MCWalkerConfiguration& w,
                               TrialWaveFunction& psi,
                               QMCHamiltonian& h,
                               RandomBase<FullPrecRealType>& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  UpdatePbyP = false;
}

SOVMCUpdateAll::~SOVMCUpdateAll() {}

void SOVMCUpdateAll::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  /* thisWalker.R will track the last accepted configuration
   * W.R will track the proposed configuration 
   *
   * upon Call:
   *   thisWalker.R,G,L,Properties must be consistent
   *   recompute flag is for pbyp driver and is not used here
   *
   * upon Return:
   *   thisWalker.R,G,L,Properties must be kept consistent
   * */
  bool updated = false;
  W.loadWalker(
      thisWalker,
      false); // W.R,G,L = thisWalker.R,G,L; false indicates W.DistTables & SK are not updated in this call. W.DistTables,SK are now stale.
  RealType logpsi_old = thisWalker.Properties(WP::LOGPSI);

  RealType invSqrtSpinMass = 1.0 / std::sqrt(spinMass);

  for (int iter = 0; iter < nSubSteps; ++iter)
  { // make a few Monte-Carlo steps to decorrelate samples without calculating observables
    makeGaussRandomWithEngine(deltaR, RandomGen); // fill deltaR
    makeGaussRandomWithEngine(deltaS, RandomGen);
    updated = false;
    if (UseDrift)
    {
      APP_ABORT("All-electron spin moves with drift are not implemented\n");
    }
    else
    {
      if (W.makeMoveAllParticles(thisWalker, deltaR, SqrtTauOverMass))
      {
        //Assume that since the spatial move is fine, the spin move will be.
        //Compute that now.
        for (int iat = 0; iat < deltaR.size(); ++iat)
          W.spins[iat] = thisWalker.spins[iat] + invSqrtSpinMass * SqrtTauOverMass[iat] * deltaS[iat];

        // W.R += dR*dt; W.DistTables,SK are updated; W.G,L are now stale
        RealType logpsi = Psi.evaluateLog(W); // update W.G,L at changed W.R; update Psi.log_real_,PhaseValue
        RealType g      = std::exp(2.0 * (logpsi - logpsi_old));
        if (RandomGen() <= g)
        {                        // move is accepted; logpsi_old and R_old are stale
          logpsi_old   = logpsi; // update logpsi_old
          thisWalker.R = W.R;    // update R_old; side effect: thisWalker.G,L,DistTables,SK,Properties are now stale
          thisWalker.spins =
              W.spins; // update R_old; side effect: thisWalker.G,L,DistTables,SK,Properties are now stale
          ++nAccept;
          updated = true;
        }
      }
    }

    if (!updated)
    { // W.R,G,L,DistTables,SK are not updated. i.e. they are still consistent
      thisWalker.Age++;
      ++nReject;
    }
  }

  if (!updated)
  {                                  // W.G and W.L have to be computed because the last move was rejected
    W.loadWalker(thisWalker, false); // move W back to last accepted configuration
    W.update();                      // update W.DistTables and SK
    logpsi_old = Psi.evaluateLog(W); // update W.G,L
  }                                  // W and logpsi_old are up-to-date at this point

  RealType eloc = H.evaluate(
      W); // calculate local energy; W.SK must be up-to-date if Coulomb interaction is used with periodic boundary. W.SK is used to calculate the long-range part of the Coulomb potential.
  W.saveWalker(thisWalker);
  thisWalker.resetProperty(logpsi_old, Psi.getPhase(),
                           eloc);               // update thisWalker::Properties[WP::LOGPSI,WP::SIGN,WP::LOCALENERGY]
  H.auxHevaluate(W, thisWalker);                // update auxiliary observables, i.e. fill H::Observables
  H.saveProperty(thisWalker.getPropertyBase()); // copy H::Observables to thisWalker::Properties
}

} // namespace qmcplusplus
