//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: DMC.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <functional>
#include <cassert>
#include <cmath>

#include "QMCDrivers/DMC/DMCBatched.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Concurrency/Info.hpp"
#include "Utilities/RunTimeManager.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
using std::placeholders::_1;
using WP = WalkerProperties::Indexes;

// clang-format off
/** Constructor maintains proper ownership of input parameters
 *
 *  Note you must call the Base constructor before the derived class sets QMCType
 */
DMCBatched::DMCBatched(QMCDriverInput&& qmcdriver_input,
                       DMCDriverInput&& input,
                       MCPopulation& pop,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       WaveFunctionPool& wf_pool,
                       Communicate* comm)
    : QMCDriverNew(std::move(qmcdriver_input), pop, psi, h, wf_pool,
                   "DMCBatched::", comm,
                   std::bind(&DMCBatched::setNonLocalMoveHandler, this, _1)),
      dmcdriver_input_(input),
      dmc_timers_("DMCBatched::")
{
  QMCType = "DMCBatched";
}
// clang-format on

void DMCBatched::setNonLocalMoveHandler(QMCHamiltonian& golden_hamiltonian)
{
  golden_hamiltonian.setNonLocalMoves(dmcdriver_input_.get_non_local_move(), qmcdriver_input_.get_tau(),
                                      dmcdriver_input_.get_alpha(), dmcdriver_input_.get_gamma());
}

void DMCBatched::resetUpdateEngines()
{
  ReportEngine PRE("DMC", "resetUpdateEngines");
  Timer init_timer;
  // Here DMC loads "Ensemble of cloned MCWalkerConfigurations"
  // I'd like to do away with this method in DMCBatched.

  // false indicates we do not support kill at node crossings.
  branch_engine_->initWalkerController(population_, dmcdriver_input_.get_reconfiguration(), false);

  estimator_manager_->reset();

  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population_.get_walkers()));

  branch_engine_->checkParameters(population_.get_num_global_walkers(), walkers);

  std::ostringstream o;
  if (dmcdriver_input_.get_reconfiguration())
    o << "  Fixed population using reconfiguration method\n";
  else
    o << "  Fluctuating population\n";
  o << "  Persistent walkers are killed after " << dmcdriver_input_.get_max_age() << " MC sweeps\n";
  o << "  BranchInterval = " << dmcdriver_input_.get_branch_interval() << "\n";
  o << "  Steps per block = " << qmcdriver_input_.get_max_steps() << "\n";
  o << "  Number of blocks = " << qmcdriver_input_.get_max_blocks() << "\n";
  app_log() << o.str() << std::endl;

  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs" << std::endl;
}

void DMCBatched::advanceWalkers(const StateForThread& sft,
                                Crowd& crowd,
                                DriverTimers& timers,
                                DMCTimers& dmc_timers,
                                ContextForSteps& step_context,
                                bool recompute)
{
  timers.buffer_timer.start();
  // We copy positions from the walkers to elec particles sets for all the crowds walkers
  // we might have received a few updates over the wire (from MPI)
  // ** None of the per walker objects have referential integrity step to step **
  crowd.loadWalkers();

  int nnode_crossing(0);
  auto& walker_twfs  = crowd.get_walker_twfs();
  auto& walkers      = crowd.get_walkers();
  auto& walker_elecs = crowd.get_walker_elecs();

  // Note this resets the identities of all the walker TWFs
  auto copyTWFFromBuffer = [](TrialWaveFunction& twf, ParticleSet& pset, MCPWalker& walker) {
    twf.copyFromBuffer(pset, walker.DataSet);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    copyTWFFromBuffer(walker_twfs[iw], walker_elecs[iw], walkers[iw]);

  timers.buffer_timer.stop();

  timers.movepbyp_timer.start();
  const int num_walkers = crowd.size();
  //This generates an entire steps worth of deltas.
  step_context.nextDeltaRs(num_walkers * sft.population.get_num_particles());
  auto it_delta_r = step_context.deltaRsBegin();

  std::vector<TrialWaveFunction::GradType> grads_now(num_walkers, TrialWaveFunction::GradType(0.0));
  std::vector<TrialWaveFunction::GradType> grads_new(num_walkers, TrialWaveFunction::GradType(0.0));
  std::vector<TrialWaveFunction::PsiValueType> ratios(num_walkers, TrialWaveFunction::PsiValueType(0.0));
  std::vector<PosType> drifts(num_walkers, 0.0);
  std::vector<RealType> log_gf(num_walkers, 0.0);
  std::vector<RealType> log_gb(num_walkers, 0.0);
  std::vector<RealType> prob(num_walkers, 0.0);

  // local list to handle accept/reject
  std::vector<bool> isAccepted;
  std::vector<std::reference_wrapper<ParticleSet>> elec_accept_list, elec_reject_list;
  isAccepted.reserve(num_walkers);
  elec_accept_list.reserve(num_walkers);
  elec_reject_list.reserve(num_walkers);

  //copy the old energies
  std::vector<FullPrecRealType> old_walker_energies(num_walkers);
  auto readOldEnergies = [](MCPWalker& walker, FullPrecRealType& old_walker_energy) {
    old_walker_energy = walker.Properties(WP::LOCALENERGY);
  };
  for (int iw = 0; iw < num_walkers; ++iw)
    readOldEnergies(walkers[iw], old_walker_energies[iw]);
  std::vector<FullPrecRealType> new_walker_energies{old_walker_energies};

  std::vector<RealType> gf_acc(num_walkers, 1.0);
  std::vector<RealType> rr_proposed(num_walkers, 0.0);
  std::vector<RealType> rr_accepted(num_walkers, 0.0);

  //  dmc_timers_.dmc_movePbyP.start();

  std::vector<int> did_walker_move(num_walkers, 0);

  for (int ig = 0; ig < step_context.get_num_groups(); ++ig)
  {
    RealType tauovermass = sft.qmcdrv_input.get_tau() * sft.population.get_ptclgrp_inv_mass()[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
    int start_index      = step_context.getPtclGroupStart(ig);
    int end_index        = step_context.getPtclGroupEnd(ig);
    for (int iat = start_index; iat < end_index; ++iat)
    {
      auto delta_r_start = it_delta_r + iat * num_walkers;
      auto delta_r_end   = delta_r_start + num_walkers;

      //This is very useful thing to be able to look at in the debugger
#ifndef NDEBUG
      std::vector<int> walkers_who_have_been_on_wire(num_walkers, 0);
      ;
      for (int iw = 0; iw < walkers.size(); ++iw)
      {
        walkers[iw].get().get_has_been_on_wire() ? walkers_who_have_been_on_wire[iw] = 1
                                                 : walkers_who_have_been_on_wire[iw] = 0;
      }
#endif
      //get the displacement
      TrialWaveFunction::flex_evalGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, grads_now);
      sft.drift_modifier.getDrifts(tauovermass, grads_now, drifts);

      std::transform(drifts.begin(), drifts.end(), delta_r_start, drifts.begin(),
                     [sqrttau](PosType& drift, PosType& delta_r) { return drift + (sqrttau * delta_r); });

      // only DMC does this
      // TODO: rr needs a real name
      std::vector<RealType> rr(num_walkers, 0.0);
      assert(rr.size() == delta_r_end - delta_r_start);
      std::transform(delta_r_start, delta_r_end, rr.begin(),
                     [tauovermass](auto& delta_r) { return tauovermass * dot(delta_r, delta_r); });

// in DMC this was done here, changed to match VMCBatched pending factoring to common source
// if (rr > m_r2max)
// {
//   ++nRejectTemp;
//   continue;
// }
#ifndef NDEBUG
      for (int i = 0; i < rr.size(); ++i)
        assert(std::isfinite(rr[i]));
#endif
      auto elecs = crowd.get_walker_elecs();
      ParticleSet::flex_makeMove(crowd.get_walker_elecs(), iat, drifts);

      TrialWaveFunction::flex_calcRatioGrad(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, ratios, grads_new);

      // This lambda is not nested thread safe due to the nreject, nnode_crossing updates
      auto checkPhaseChanged = [&sft, &iat, &crowd, &nnode_crossing](TrialWaveFunction& twf, ParticleSet& elec,
                                                                     int& is_reject) {
        if (sft.branch_engine.phaseChanged(twf.getPhaseDiff()))
        {
          crowd.incReject();
          ++nnode_crossing;
          elec.rejectMove(iat);
          twf.rejectMove(iat);
          is_reject = 1;
        }
        else
          is_reject = 0;
      };

      // Hopefully a phase change doesn't make any of these transformations fail.
      std::vector<int> rejects(num_walkers); // instead of std::vector<bool>
      for (int iw = 0; iw < num_walkers; ++iw)
      {
        checkPhaseChanged(walker_twfs[iw], walker_elecs[iw], rejects[iw]);
        //This is just convenient to do here
        rr_proposed[iw] += rr[iw];
      }

      std::transform(delta_r_start, delta_r_end, log_gf.begin(), [](auto& delta_r) {
        constexpr RealType mhalf(-0.5);
        return mhalf * dot(delta_r, delta_r);
      });

      sft.drift_modifier.getDrifts(tauovermass, grads_new, drifts);

      std::transform(crowd.beginElectrons(), crowd.endElectrons(), drifts.begin(), drifts.begin(),
                     [iat](auto& elecs, auto& drift) { return elecs.get().R[iat] - elecs.get().activePos - drift; });

      std::transform(drifts.begin(), drifts.end(), log_gb.begin(),
                     [oneover2tau](auto& drift) { return -oneover2tau * dot(drift, drift); });

      for (int iw = 0; iw < num_walkers; ++iw)
        prob[iw] = std::norm(ratios[iw]) * std::exp(log_gb[iw] - log_gf[iw]);

      isAccepted.clear();
      elec_accept_list.clear();
      elec_reject_list.clear();

      for (int iw = 0; iw < num_walkers; ++iw)
      {
        // Allows the rng to be scrutinized easily in debugger
        auto the_rng = step_context.get_random_gen()();
        if ((!rejects[iw]) && prob[iw] >= std::numeric_limits<RealType>::epsilon() &&
            the_rng < prob[iw])
        {
          did_walker_move[iw] += 1;
          crowd.incAccept();
          isAccepted.push_back(true);
          elec_accept_list.push_back(crowd.get_walker_elecs()[iw]);
          rr_accepted[iw] += rr[iw];
          gf_acc[iw] *= prob[iw];
        }
        else
        {
          crowd.incReject();
          isAccepted.push_back(false);
          elec_reject_list.push_back(crowd.get_walker_elecs()[iw]);
        }
      }

      TrialWaveFunction::flex_accept_rejectMove(crowd.get_walker_twfs(), crowd.get_walker_elecs(), iat, isAccepted,
                                                true);

      ParticleSet::flex_acceptMove(elec_accept_list, iat, true);
      ParticleSet::flex_rejectMove(elec_reject_list, iat);
    }
  }

  TrialWaveFunction::flex_completeUpdates(crowd.get_walker_twfs());
  ParticleSet::flex_donePbyP(crowd.get_walker_elecs());
  //dmc_timers.dmc_movePbyP.stop();
  timers.movepbyp_timer.stop();

  //To use the flex interfaces we have to build RefVectors for walker that moved and walkers that didn't

  auto& walker_hamiltonians  = crowd.get_walker_hamiltonians();
  auto& walker_mcp_wfbuffers = crowd.get_mcp_wfbuffers();

  DMCPerWalkerRefRefs per_walker_ref_refs{walkers,
                                          walker_twfs,
                                          walker_hamiltonians,
                                          walker_elecs,
                                          walker_mcp_wfbuffers,
                                          old_walker_energies,
                                          new_walker_energies,
                                          rr_proposed,
                                          rr_accepted,
                                          gf_acc};

  MovedStalled these = buildMovedStalled(did_walker_move, per_walker_ref_refs);

  handleMovedWalkers(these.moved, sft, timers);
  handleStalledWalkers(these.stalled, sft);

  try
  {
    dmc_timers.tmove_timer.start();
    std::vector<int> walker_non_local_moves_accepted(
        QMCHamiltonian::flex_makeNonLocalMoves(crowd.get_walker_hamiltonians(), crowd.get_walker_elecs()));

    //could be premature optimization
    int num_moved_nonlocal   = 0;
    int total_moved_nonlocal = 0;
    for (int iw = 0; iw < walkers.size(); ++iw)
    {
      if (walker_non_local_moves_accepted[iw] > 0)
      {
        num_moved_nonlocal++;
        total_moved_nonlocal += walker_non_local_moves_accepted[iw];
        crowd.incNonlocalAccept();
      }
    }

    if (num_moved_nonlocal > 0)
    {
      DMCPerWalkerRefs moved_nonlocal(num_moved_nonlocal);

      for (int iw = 0; iw < these.moved.walkers.size(); ++iw)
      {
        if (walker_non_local_moves_accepted[iw] > 0)
        {
          moved_nonlocal.walkers.push_back(walkers[iw]);
          moved_nonlocal.walker_twfs.push_back(walker_twfs[iw]);
          moved_nonlocal.walker_elecs.push_back(walker_elecs[iw]);
          moved_nonlocal.walker_hamiltonians.push_back(walker_hamiltonians[iw]);
          moved_nonlocal.walker_mcp_wfbuffers.push_back(walker_mcp_wfbuffers[iw]);
        }
      }
      TrialWaveFunction::flex_updateBuffer(moved_nonlocal.walker_twfs, moved_nonlocal.walker_elecs,
                                           moved_nonlocal.walker_mcp_wfbuffers);
      ParticleSet::flex_saveWalker(moved_nonlocal.walker_elecs, moved_nonlocal.walkers);
    }
  }
  catch (const std::out_of_range& exc)
  {
    std::cout << "Out of range error in non local move updates: " << exc.what() << std::endl;
  }
  catch (const std::exception& exc)
  {
    std::cout << "Serious issue in flex_make_NonLocalMoves(): " << exc.what() << std::endl;
  }

  dmc_timers.tmove_timer.stop();

  setMultiplicities(sft.dmcdrv_input, walkers, step_context.get_random_gen());
}

DMCBatched::MovedStalled DMCBatched::buildMovedStalled(const std::vector<int>& did_walker_move,
                                                       const DMCPerWalkerRefRefs& refs)
{
  int num_walkers = refs.walkers.size();
  int num_moved   = 0;
  for (int iw = 0; iw < num_walkers; ++iw)
  {
    if (did_walker_move[iw] > 0)
      num_moved++;
  }

  MovedStalled these(num_walkers, num_moved);
  for (int iw = 0; iw < num_walkers; ++iw)
  {
    if (did_walker_move[iw] > 0)
    {
      assert(refs.rr_accepted[iw] > 0);
      these.moved.walkers.push_back(refs.walkers[iw]);
      these.moved.walker_twfs.push_back(refs.walker_twfs[iw]);
      these.moved.walker_hamiltonians.push_back(refs.walker_hamiltonians[iw]);
      these.moved.walker_elecs.push_back(refs.walker_elecs[iw]);
      these.moved.walker_mcp_wfbuffers.push_back(refs.walker_mcp_wfbuffers[iw]);
      these.moved.old_energies.push_back(refs.old_energies[iw]);
      these.moved.new_energies.push_back(refs.new_energies[iw]);
      these.moved.rr_proposed.push_back(refs.rr_proposed[iw]);
      these.moved.rr_accepted.push_back(refs.rr_accepted[iw]);
      these.moved.gf_accs.push_back(refs.gf_accs[iw]);
    }
    else
    {
      assert(refs.rr_accepted[iw] == 0.0);
      these.stalled.walkers.push_back(refs.walkers[iw]);
      these.stalled.walker_twfs.push_back(refs.walker_twfs[iw]);
      these.stalled.walker_hamiltonians.push_back(refs.walker_hamiltonians[iw]);
      these.stalled.walker_elecs.push_back(refs.walker_elecs[iw]);
      these.stalled.walker_mcp_wfbuffers.push_back(refs.walker_mcp_wfbuffers[iw]);
      these.stalled.old_energies.push_back(refs.old_energies[iw]);
      these.stalled.new_energies.push_back(refs.new_energies[iw]);
      these.stalled.rr_proposed.push_back(refs.rr_proposed[iw]);
      these.stalled.rr_accepted.push_back(refs.rr_accepted[iw]);
      these.stalled.gf_accs.push_back(refs.gf_accs[iw]);
    }
  }
  return these;
}

void DMCBatched::handleMovedWalkers(DMCPerWalkerRefs& moved, const StateForThread& sft, DriverTimers& timers)
{
  if (moved.walkers.size() > 0)
  {
    timers.buffer_timer.start();
    TrialWaveFunction::flex_updateBuffer(moved.walker_twfs, moved.walker_elecs, moved.walker_mcp_wfbuffers);
    std::for_each(moved.walkers.begin(), moved.walkers.end(), [](MCPWalker& walker) { walker.Age = 0; });
    ParticleSet::flex_saveWalker(moved.walker_elecs, moved.walkers);
    timers.buffer_timer.stop();
    timers.hamiltonian_timer.start();
    // std::vector<QMCHamiltonian::FullPrecRealType> local_energies(QMCHamiltonian::flex_evaluate(moved.walker_hamiltonians, moved.walker_elecs));
    std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
        QMCHamiltonian::flex_evaluateWithToperator(moved.walker_hamiltonians, moved.walker_elecs));
    timers.hamiltonian_timer.stop();

    auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto local_energy, auto rr_acc,
                                   auto rr_prop) {
      walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy, rr_acc, rr_prop, 1.0);
    };

    for (int iw = 0; iw < moved.walkers.size(); ++iw)
    {
      assert(moved.rr_proposed[iw] > 0);
      resetSigNLocalEnergy(moved.walkers[iw], moved.walker_twfs[iw], local_energies[iw], moved.rr_accepted[iw],
                           moved.rr_proposed[iw]);
      // this might mean new_energies are actually unneeded which would be nice.
      moved.new_energies[iw]         = local_energies[iw];
      FullPrecRealType branch_weight = sft.branch_engine.branchWeight(moved.new_energies[iw], moved.old_energies[iw]);
      moved.walkers[iw].get().Weight *= branch_weight;
    }
    timers.collectables_timer.start();
    auto evaluateNonPhysicalHamiltonianElements = [](QMCHamiltonian& ham, ParticleSet& pset, MCPWalker& walker) {
      ham.auxHevaluate(pset, walker);
    };
    for (int iw = 0; iw < moved.walkers.size(); ++iw)
      evaluateNonPhysicalHamiltonianElements(moved.walker_hamiltonians[iw], moved.walker_elecs[iw], moved.walkers[iw]);

    auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
      ham.saveProperty(walker.getPropertyBase());
    };
    for (int iw = 0; iw < moved.walkers.size(); ++iw)
      savePropertiesIntoWalker(moved.walker_hamiltonians[iw], moved.walkers[iw]);
    timers.collectables_timer.stop();
  }
}

void DMCBatched::handleStalledWalkers(DMCPerWalkerRefs& stalled, const StateForThread& sft)
{
  for (int iw = 0; iw < stalled.walkers.size(); ++iw)
  {
    std::cout << "A walker has stalled.\n";
    TrialWaveFunction::flex_updateBuffer(stalled.walker_twfs, stalled.walker_elecs, stalled.walker_mcp_wfbuffers);
    std::for_each(stalled.walkers.begin(), stalled.walkers.end(), [](MCPWalker& walker) { walker.Age = 0; });
    ParticleSet::flex_saveWalker(stalled.walker_elecs, stalled.walkers);

    MCPWalker& stalled_walker = stalled.walkers[iw];
    stalled_walker.Age++;
    stalled_walker.Properties(WP::R2ACCEPTED) = 0.0;
    FullPrecRealType wtmp                     = stalled_walker.Weight;
    // TODO: fix this walker.Weight twiddle for rejectedMove
    stalled_walker.Weight                      = 0.0;
    QMCHamiltonian& stalled_walker_hamiltonian = stalled.walker_hamiltonians[iw];
    ParticleSet& stalled_particle_set          = stalled.walker_elecs[iw];
    stalled_walker_hamiltonian.rejectedMove(stalled_particle_set, stalled_walker);
    stalled_walker.Weight                       = wtmp;
    FullPrecRealType& stalled_new_walker_energy = stalled.new_energies[iw];
    FullPrecRealType& stalled_old_walker_energy = stalled.old_energies[iw];
    stalled_new_walker_energy                   = stalled_old_walker_energy;
    stalled.gf_accs[iw].get()                   = 1.0;
    stalled_walker.Weight *= sft.branch_engine.branchWeight(stalled_new_walker_energy, stalled_old_walker_energy);
  }
}

void DMCBatched::setMultiplicities(const DMCDriverInput& dmcdriver_input,
                                   RefVector<MCPWalker>& walkers,
                                   RandomGenerator_t& rng)
{
  auto setMultiplicity = [&dmcdriver_input, &rng](MCPWalker& walker) {
    constexpr FullPrecRealType onehalf(0.5);
    constexpr FullPrecRealType cone(1);
    walker.Multiplicity = walker.Weight;
    if (walker.Age > dmcdriver_input.get_max_age())
      walker.Multiplicity = std::min(onehalf, walker.Weight);
    else if (walker.Age > 0)
      walker.Multiplicity = std::min(cone, walker.Weight);
    walker.Multiplicity += rng();
  };

  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    setMultiplicity(walkers[iw]);
  }
}

void DMCBatched::runDMCStep(int crowd_id,
                            const StateForThread& sft,
                            DriverTimers& timers,
                            DMCTimers& dmc_timers,
                            UPtrVector<ContextForSteps>& context_for_steps,
                            UPtrVector<Crowd>& crowds)
{
  Crowd& crowd = *(crowds[crowd_id]);
  if (crowd.size() == 0)
    return;
  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());

  int max_steps = sft.qmcdrv_input.get_max_steps();
  // This is migraine inducing here and in the original driver, I believe they are the same in
  // VMC(Batched)/DMC(Batched) needs another check and unit test
  bool is_recompute_block =
      sft.recomputing_blocks ? (1 + sft.block) % sft.qmcdrv_input.get_blocks_between_recompute() == 0 : false;
  IndexType step           = sft.step;
  bool recompute_this_step = (is_recompute_block && (step + 1) == max_steps);
  advanceWalkers(sft, crowd, timers, dmc_timers, *context_for_steps[crowd_id], recompute_this_step);
}

void DMCBatched::process(xmlNodePtr node)
{
  QMCDriverNew::AdjustedWalkerCounts awc =
      adjustGlobalWalkerCount(myComm->size(), myComm->rank(), qmcdriver_input_.get_total_walkers(),
                              qmcdriver_input_.get_walkers_per_rank(), dmcdriver_input_.get_reserve(),
                              qmcdriver_input_.get_num_crowds());

  Base::startup(node, awc);
}

bool DMCBatched::run()
{
  resetUpdateEngines();
  IndexType num_blocks = qmcdriver_input_.get_max_blocks();

  estimator_manager_->start(num_blocks);
  StateForThread dmc_state(qmcdriver_input_, dmcdriver_input_, *drift_modifier_, *branch_engine_, population_);

  LoopTimer<> dmc_loop;

  int sample = 0;
  RunTimeControl<> runtimeControl(run_time_manager, MaxCPUSecs);

  { // walker initialization
    ScopedTimer local_timer(&(timers_.init_walkers_timer));
    TasksOneToOne<> section_start_task(crowds_.size());
    section_start_task(initialLogEvaluation, std::ref(crowds_), std::ref(step_contexts_));
  }

  TasksOneToOne<> crowd_task(crowds_.size());

  for (int block = 0; block < num_blocks; ++block)
  {
    dmc_loop.start();
    estimator_manager_->startBlock(qmcdriver_input_.get_max_steps());

    dmc_state.recalculate_properties_period = (qmc_driver_mode_[QMC_UPDATE_MODE])
        ? qmcdriver_input_.get_recalculate_properties_period()
        : (qmcdriver_input_.get_max_blocks() + 1) * qmcdriver_input_.get_max_steps();
    dmc_state.recomputing_blocks            = qmcdriver_input_.get_blocks_between_recompute();

    for (UPtr<Crowd>& crowd : crowds_)
      crowd->startBlock(qmcdriver_input_.get_max_steps());

    for (int step = 0; step < qmcdriver_input_.get_max_steps(); ++step)
    {
      ScopedTimer local_timer(&(timers_.run_steps_timer));
      dmc_state.step = step;
      crowd_task(runDMCStep, dmc_state, timers_, dmc_timers_, std::ref(step_contexts_), std::ref(crowds_));

      // Accumulate on the whole population
      // But it is now visible in the algorithm not hidden in the BranchEngine::branch.
      // \todo make task block
      // probably something smart can be done to include reduction over crowds below
      for (UPtr<Crowd>& crowd_ptr : crowds_)
      {
        Crowd& crowd_ref = *crowd_ptr;
        if (crowd_ref.size() > 0)
          crowd_ref.accumulate(population_.get_num_global_walkers());
      }
      
      branch_engine_->branch(step, population_);

      for (UPtr<Crowd>& crowd_ptr : crowds_)
        crowd_ptr->clearWalkers();

      population_.distributeWalkers(crowds_);
    }
    endBlock();
  }
  return finalize(num_blocks, true);
}

} // namespace qmcplusplus
