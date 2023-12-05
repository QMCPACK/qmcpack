//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "RMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Concurrency/OpenMP.h"
#include "Configuration.h"
#include "Particle/Reptile.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////////
//  This driver is strongly based off the method used by Lucas Wagner in QWalk.
//
//  This driver proposes an "all-electron" configuration by performing N single particle moves,
//  accepting or rejecting at each step.  Umrigar's scaled drift is used to generate the VMC configurations.
//
//  The configuration is then accepted/rejected based only on the symmetrized DMC action, which
//  amounts to Maroni/Baroni's method.  Note that energy filtering (like in DMC) is used here as well.
//
//  Until serious discrepencies are detected, it is strongly advised that this driver be used over the
//  RMCUpdateAll, as the time step errors seem to be really reduced using this method.
//////////////////////////////////////////////////////////////////////////////


namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/// Constructor.
RMCUpdatePbyPWithDrift::RMCUpdatePbyPWithDrift(MCWalkerConfiguration& w,
                                               TrialWaveFunction& psi,
                                               QMCHamiltonian& h,
                                               RandomBase<FullPrecRealType>& rg,
                                               std::vector<int> act,
                                               std::vector<int> tp)
    : QMCUpdateBase(w, psi, h, rg),
      Action(act),
      TransProb(tp),
      advance_timer_(createGlobalTimer("RMCUpdatePbyP::advance", timer_level_medium)),
      movepbyp_timer_(createGlobalTimer("RMCUpdatePbyP::movePbyP", timer_level_medium)),
      update_mbo_timer_(createGlobalTimer("RMCUpdatePbyP::updateMBO", timer_level_medium)),
      energy_timer_(createGlobalTimer("RMCUpdatePbyP::energy", timer_level_medium))
{
  scaleDrift = false;
  actionType = SYM_ACTION;
}

RMCUpdatePbyPWithDrift::~RMCUpdatePbyPWithDrift() {}


void RMCUpdatePbyPWithDrift::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP = true;

  for (; it != it_end; ++it)
  {
    Walker_t& awalker = **it; //W.reptile->getHead();
    W.R               = awalker.R;
    W.update(true);
    W.donePbyP();
    if (awalker.DataSet.size())
      awalker.DataSet.clear();
    awalker.DataSet.rewind();
    Psi.registerData(W, awalker.DataSet);
    awalker.DataSet.allocate();
    Psi.copyFromBuffer(W, awalker.DataSet);
    Psi.evaluateLog(W);
    RealType logpsi = Psi.updateBuffer(W, awalker.DataSet, false);
    awalker.G       = W.G;
    awalker.L       = W.L;
    RealType eloc   = H.evaluate(W);
    awalker.resetProperty(logpsi, Psi.getPhase(), eloc);
  }
}

bool RMCUpdatePbyPWithDrift::put(xmlNodePtr cur)
{
  QMCUpdateBase::put(cur);

  ParameterSet m_param;
  bool usedrift      = true;
  std::string action = "SLA";
  m_param.add(usedrift, "useDrift");
  m_param.add(action, "Action");
  m_param.add(equilSteps, "equilsteps");
  m_param.add(equilSteps, "equilSteps");
  m_param.put(cur);

  if (usedrift == true)
  {
    if (omp_get_thread_num() == 0)
      app_log() << "  Using Umrigar scaled drift\n";
  }
  else
  {
    if (omp_get_thread_num() == 0)
      app_log() << "  Using non-scaled drift\n";
  }

  if (action == "DMC")
  {
    actionType = DMC_ACTION;
    if (omp_get_thread_num() == 0)
      app_log() << "  Using DMC link-action\n";
  }
  else
  {
    if (omp_get_thread_num() == 0)
      app_log() << "  Using Symmetrized Link-Action\n";
  }

  return true;
}
void RMCUpdatePbyPWithDrift::advanceWalkersVMC()
{
  advance_timer_.start();
  Walker_t& curhead = W.reptile->getHead();
  Walker_t prophead(curhead);
  Walker_t::WFBuffer_t& w_buffer(prophead.DataSet);
  W.loadWalker(prophead, true);
  Psi.copyFromBuffer(W, w_buffer);
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  int nAcceptTemp(0);
  //copy the old energy and scale factor of drift
  RealType eold(prophead.Properties(WP::LOCALENERGY));
  RealType vqold(prophead.Properties(WP::DRIFTSCALE));
  RealType enew(eold);
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  movepbyp_timer_.start();
  for (int ig = 0; ig < W.groups(); ++ig) //loop over species
  {
    RealType tauovermass = Tau * MassInvS[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
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
        W.accept_rejectMove(iat, false);
        continue;
      }
      ValueType ratio = Psi.calcRatioGrad(W, iat, grad_iat);
      //node is crossed reject the move
      if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nNodeCrossing;
        W.accept_rejectMove(iat, false);
        Psi.rejectMove(iat);
      }
      else
      {
        RealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
        //Use the force of the particle iat
        DriftModifier->getDrift(tauovermass, grad_iat, dr);
        dr               = W.R[iat] - W.getActivePos() - dr;
        RealType logGb   = -oneover2tau * dot(dr, dr);
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
          Psi.rejectMove(iat);
        }
        W.accept_rejectMove(iat, is_accepted);
      }
    }
  }
  movepbyp_timer_.stop();
  Psi.completeUpdates();
  W.donePbyP();

  if (nAcceptTemp > 0)
  {
    //need to overwrite the walker properties
    MCWalkerConfiguration::Walker_t& newhead(W.reptile->getNewHead());
    update_mbo_timer_.start();
    prophead.Age    = 0;
    prophead.R      = W.R;
    RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
    W.saveWalker(prophead);
    update_mbo_timer_.stop();
    energy_timer_.start();
    enew = H.evaluate(W);
    energy_timer_.stop();
    prophead.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 0.0);
    prophead.Weight = 1.0;
    H.auxHevaluate(W, prophead, true, false); //evaluate properties but not collectables.
    H.saveProperty(prophead.getPropertyBase());
    newhead = prophead;
    nAccept++;
  }
  else
  {
    //all moves are rejected: does not happen normally with reasonable wavefunctions
    curhead.Age++;
    curhead.Properties(WP::R2ACCEPTED) = 0.0;
    //weight is set to 0 for traces
    // consistent w/ no evaluate/auxHevaluate
    RealType wtmp  = prophead.Weight;
    curhead.Weight = 0.0;
    H.rejectedMove(W, curhead);
    curhead.Weight = wtmp;
    ++nAllRejected;
    nReject++;
  }
  Walker_t& centerbead = W.reptile->getCenter();
  W.loadWalker(centerbead, true);
  W.update();
  H.auxHevaluate(W, centerbead); //evaluate collectables but not properties.
  // Traces->buffer_sample();
}


void RMCUpdatePbyPWithDrift::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  IndexType initsteps = W.reptile->nbeads * 2;

  vmcSteps = W.reptile->nbeads + 1;

  for (int n = 0; n < initsteps; n++)
    advanceWalkersVMC();
}

void RMCUpdatePbyPWithDrift::advanceWalkersRMC()
{
  Walker_t& curhead = W.reptile->getHead();
  Walker_t prophead(curhead);
  Walker_t::WFBuffer_t& w_buffer(prophead.DataSet);
  W.loadWalker(prophead, true);
  Psi.copyFromBuffer(W, w_buffer);

  makeGaussRandomWithEngine(deltaR, RandomGen);
  int nAcceptTemp(0);
  //copy the old energy and scale factor of drift
  RealType eold(prophead.Properties(WP::LOCALENERGY));
  RealType vqold(prophead.Properties(WP::DRIFTSCALE));
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  movepbyp_timer_.start();
  for (int ig = 0; ig < W.groups(); ++ig) //loop over species
  {
    RealType tauovermass = Tau * MassInvS[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
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
        W.accept_rejectMove(iat, false);
        continue;
      }
      ValueType ratio = Psi.calcRatioGrad(W, iat, grad_iat);
      //node is crossed reject the move
      if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nNodeCrossing;
        W.accept_rejectMove(iat, false);
        Psi.rejectMove(iat);
      }
      else
      {
        RealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
        //Use the force of the particle iat
        DriftModifier->getDrift(tauovermass, grad_iat, dr);
        dr               = W.R[iat] - W.getActivePos() - dr;
        RealType logGb   = -oneover2tau * dot(dr, dr);
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
          Psi.rejectMove(iat);
        }
        W.accept_rejectMove(iat, is_accepted);
      }
    }
  }
  movepbyp_timer_.stop();
  Psi.completeUpdates();
  W.donePbyP();
  // In the rare case that all proposed moves fail, we bounce.
  if (nAcceptTemp == 0)
  {
    ++nReject;
    H.rejectedMove(W, prophead);
    curhead.Age += 1;
    W.reptile->flip();
  }
  prophead.R      = W.R;
  RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
  W.saveWalker(prophead);
  Walker_t &lastbead(W.reptile->getTail()), nextlastbead(W.reptile->getNext());
  RealType eloc = H.evaluate(W);
  RealType dS   = branchEngine->DMCLinkAction(eloc, curhead.Properties(WP::LOCALENERGY)) -
      branchEngine->DMCLinkAction(lastbead.Properties(WP::LOCALENERGY), nextlastbead.Properties(WP::LOCALENERGY));
  RealType acceptProb = std::min((RealType)1.0, std::exp(-dS));
  if ((RandomGen() <= acceptProb) || (prophead.Age >= MaxAge || lastbead.Age >= MaxAge))
  {
    MCWalkerConfiguration::Walker_t& overwriteWalker(W.reptile->getNewHead());
    if (curhead.Age >= MaxAge || lastbead.Age >= MaxAge)
      app_log() << "\tForce Acceptance...\n";

    prophead.Properties(WP::LOCALENERGY) = eloc;
    prophead.Properties(WP::R2ACCEPTED)  = rr_accepted;
    prophead.Properties(WP::R2PROPOSED)  = rr_proposed;
    H.auxHevaluate(W, prophead, true, false); //evaluate properties? true.  collectables?  false.
    H.saveProperty(prophead.getPropertyBase());
    prophead.Age    = 0;
    overwriteWalker = prophead;
    ++nAccept;
  }
  else
  {
    ++nReject;
    H.rejectedMove(W, prophead);
    curhead.Properties(WP::R2ACCEPTED) = 0;
    curhead.Properties(WP::R2PROPOSED) = rr_proposed;
    curhead.Age += 1;
    W.reptile->flip();
    // return;
  }
  //Collectables should be evaluated on center bead.  Here we go.
  Walker_t& centerbead = W.reptile->getCenter();
  W.loadWalker(centerbead, true);
  W.update();                                 //Called to recompute S(k) and distance tables.
  H.auxHevaluate(W, centerbead, false, true); //evaluate properties?  false.  Collectables?  true.
}

void RMCUpdatePbyPWithDrift::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  //empty function to
}

void RMCUpdatePbyPWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool init)
{
  if (init == true)
    advanceWalkersVMC();
  else
    advanceWalkersRMC();
}

void RMCUpdatePbyPWithDrift::accumulate(WalkerIter_t it, WalkerIter_t it_end) { Estimators->accumulate(W, it, it_end); }

} // namespace qmcplusplus
