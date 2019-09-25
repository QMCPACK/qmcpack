//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/DMC/DMCBatched.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
/** Constructor maintains proper ownership of input parameters
   */
DMCBatched::DMCBatched(QMCDriverInput&& qmcdriver_input,
                       DMCDriverInput&& input,
                       MCPopulation&& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& wf_pool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), std::move(pop), psi, h, wf_pool, "DMCBatched::", comm), 
      dmcdriver_input_(input)
{
  QMCType = "DMCBatched";
}

QMCTraits::IndexType DMCBatched::calc_default_local_walkers(IndexType walkers_per_rank)
{
  checkNumCrowdsLTNumThreads();
  int num_threads(Concurrency::maxThreads<>());
  IndexType rw = walkers_per_rank; //qmcdriver_input_.get_walkers_per_rank();
  if (num_crowds_ == 0)
    num_crowds_ = std::min(num_threads, rw);
  walkers_per_crowd_      = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;

  IndexType local_walkers = walkers_per_crowd_ * num_crowds_;
  population_.set_num_local_walkers(local_walkers);
  population_.set_num_global_walkers(local_walkers * population_.get_num_ranks());
  if (rw != qmcdriver_input_.get_walkers_per_rank())
    app_warning() << "DMCBatched driver has adjusted walkers per rank to: " << local_walkers << '\n';

  app_log() << "DMCBatched walkers per crowd " << walkers_per_crowd_ << std::endl;
  return local_walkers;
}

void DMCBatched::resetUpdateEngines()
{
  ReportEngine PRE("DMC", "resetUpdateEngines");
  Timer init_timer;
  // Here DMC loads "Ensemble of cloned MCWalkerConfigurations"
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population_.get_walkers()));
  branchEngine->checkParameters(population_.get_num_global_walkers(), walkers);

  std::ostringstream o;
  if (dmcdriver_input_.get_reconfiguration())
    o << "  Fixed population using reconfiguration method\n";
  else
    o << "  Fluctuating population\n";
  o << "  Persistent walkers are killed after " << dmcdriver_input_.get_mover_max_age() << " MC sweeps\n";
  o << "  BranchInterval = " << dmcdriver_input_.get_branch_interval() << "\n";
  o << "  Steps per block = " << qmcdriver_input_.get_max_steps() << "\n";
  o << "  Number of blocks = " << qmcdriver_input_.get_max_blocks() << "\n";
  app_log() << o.str() << std::endl;

  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs" << std::endl;
}

void DMCBatched::advanceWalkers(const StateForThread& sft, Crowd& crowd, DriverTimers& timers, ContextForSteps& step_context, bool recompute)
{
  timers.buffer_timer.start();
  crowd.loadWalkers();

  // Consider favoring lambda followed by for over walkers
  // more compact, descriptive and less error prone.
  auto& walker_twfs      = crowd.get_walker_twfs();
  auto& walkers          = crowd.get_walkers();
  auto& walker_elecs     = crowd.get_walker_elecs();
  auto copyTWFFromBuffer = [](TrialWaveFunction& twf, ParticleSet& pset, MCPWalker& walker) {
    twf.copyFromBuffer(pset, walker.DataSet);
  };
  timers.buffer_timer.stop();

  timers.movepbyp_timer.start();
  int num_walkers = crowd.size();
  //This generates an entire steps worth of deltas.
  step_context.nextDeltaRs();
  auto it_delta_r = step_context.deltaRsBegin();
  std::vector<PosType> drifts(num_walkers);

  // local list to handle accept/reject
  std::vector<std::reference_wrapper<ParticleSet>> elec_accept_list, elec_reject_list;
  std::vector<std::reference_wrapper<TrialWaveFunction>> twf_accept_list, twf_reject_list;
  elec_accept_list.reserve(num_walkers);
  elec_reject_list.reserve(num_walkers);
  twf_accept_list.reserve(num_walkers);
  twf_reject_list.reserve(num_walkers);

//copy the old energy
  std::vector<FullPrecRealType> old_walker_energies(num_walkers);
  auto setOldEnergies = [](MCPWalker& walker, FullPrecRealType& old_walker_energy) {
                          old_walker_energy = walker.Properties(LOCALENERGY);
                        };
  for(int iw = 0; iw < num_walkers; ++iw)
    setOldEnergies(walkers[iw], old_walker_energies[iw]);
  std::vector<FullPrecRealType> new_walker_energies(old_walker_energies);

  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  RealType gf_acc      = 1.0;

  timers.myTimers[DMC_movePbyP]->start();
  for (int ig = 0; ig < W.groups(); ++ig) //loop over species
  {
    RealType tauovermass = Tau * MassInvS[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
    for (int iat = W.first(ig); iat < W.last(ig); ++iat)
    {
      W.setActive(iat);
      //get the displacement
      GradType grad_iat = Psi.evalGrad(W, iat);
      PosType dr;
      DriftModifier->getDrift(tauovermass, grad_iat, dr);
      dr += sqrttau * deltaR[iat];
      RealType rr = tauovermass * dot(deltaR[iat], deltaR[iat]);
      rr_proposed += rr;
      if (rr > m_r2max)
      {
        ++nRejectTemp;
        continue;
      }
      if (!W.makeMoveAndCheck(iat, dr))
        continue;
      RealType ratio = Psi.ratioGrad(W, iat, grad_iat);
      //node is crossed reject the move
      if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nRejectTemp;
        ++nNodeCrossing;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
      else
      {
        FullPrecRealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
        //Use the force of the particle iat
        DriftModifier->getDrift(tauovermass, grad_iat, dr);
        dr                     = W.R[iat] - W.activePos - dr;
        FullPrecRealType logGb = -oneover2tau * dot(dr, dr);
        RealType prob          = ratio * ratio * std::exp(logGb - logGf);
        if (RandomGen() < prob)
        {
          ++nAcceptTemp;
          Psi.acceptMove(W, iat);
          W.acceptMove(iat);
          rr_accepted += rr;
          gf_acc *= prob; //accumulate the ratio
        }
        else
        {
          ++nRejectTemp;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
  }
  Psi.completeUpdates();
  W.donePbyP();
  myTimers[DMC_movePbyP]->stop();

  if (nAcceptTemp > 0)
  {
    //need to overwrite the walker properties
    myTimers[DMC_buffer]->start();
    thisWalker.Age  = 0;
    RealType logpsi = Psi.updateBuffer(W, w_buffer, recompute);
    W.saveWalker(thisWalker);
    myTimers[DMC_buffer]->stop();
    myTimers[DMC_hamiltonian]->start();
    enew = H.evaluateWithToperator(W);
    myTimers[DMC_hamiltonian]->stop();
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 1.0);
    thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    myTimers[DMC_collectables]->start();
    H.auxHevaluate(W, thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());
    myTimers[DMC_collectables]->stop();
  }
  else
  {
    //all moves are rejected: does not happen normally with reasonable wavefunctions
    thisWalker.Age++;
    thisWalker.Properties(R2ACCEPTED) = 0.0;
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
  myTimers[DMC_tmoves]->start();
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
  myTimers[DMC_tmoves]->stop();
  nAccept += nAcceptTemp;
  nReject += nRejectTemp;

  setMultiplicity(thisWalker);

}

void DMCBatched::runDMCStep(int crowd_id,
                            const StateForThread& sft,
                            DriverTimers& timers,
                            std::vector<std::unique_ptr<ContextForSteps>>& context_for_steps,
                            std::vector<std::unique_ptr<Crowd>>& crowds)
{
  Crowd& crowd = *(crowds[crowd_id]);
  int max_steps = sft.qmcdrv_input.get_max_steps();
  // This is migraine inducing here and in the original driver, I believe they are the same in
  // VMC(Batched)/DMC(Batched) needs another check and unit test
  bool is_recompute_block =
      sft.recomputing_blocks ? (1 + sft.block) % sft.qmcdrv_input.get_blocks_between_recompute() == 0 : false;
  IndexType step = sft.step;
  bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps);
  advanceWalkers(sft, crowd, timers, *context_for_steps[crowd_id], recompute_this_step);
}

bool DMCBatched::run()
{
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();

  estimator_manager_->setCollectionMode(true);
  estimator_manager_->start(num_blocks);
  StateForThread dmc_state(qmcdriver_input_, dmcdriver_input_, *drift_modifier_, population_);

  LoopTimer dmc_loop;

  int sample             = 0;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);

  TasksOneToOne<> crowd_task(num_crowds_);

  for (int block = 0; block < num_blocks; ++block)
  {
    dmc_loop.start();
    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    dmc_state.recalculate_properties_period = (qmc_driver_mode_[QMC_UPDATE_MODE]) ? qmcdriver_input_.get_recalculate_properties_period() : (qmcdriver_input_.get_max_blocks() + 1) * qmcdriver_input_.get_max_steps();
    dmc_state.recomputing_blocks = qmcdriver_input_.get_blocks_between_recompute();

    for (auto& crowd : crowds_)
      crowd->startBlock(qmcdriver_input_.get_max_steps());

    for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      ScopedTimer local_timer(&(timers_.run_steps_timer));
      dmc_state.step = step;
      crowd_task(runDMCStep, dmc_state, timers_, std::ref(step_contexts_), std::ref(crowds_));
    }
  }
  return false;
}

} // namespace qmcplusplus
