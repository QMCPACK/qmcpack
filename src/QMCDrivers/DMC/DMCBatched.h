//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from DMC.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMCBATCHED_H
#define QMCPLUSPLUS_DMCBATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/DMC/DMCDriverInput.h"
#include "QMCDrivers/MCPopulation.h"
#include "Particle/MCCoords.hpp"
#include "WalkerLogManager.h"
#include "RunTimeManager.h"

namespace qmcplusplus
{
class DriverModifierBase;
class WalkerControl;
class SFNBranch;

namespace testing
{
class DMCBatchedTest;
class DMCBatchedTestAccessor;
} // namespace testing

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a DMC using particle-by-particle threaded and batched moves.
 */
class DMCBatched : public QMCDriverNew
{
public:
  using FullPrecRealType  = QMCTraits::FullPrecRealType;
  using PosType           = QMCTraits::PosType;
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos;
  /** To avoid 10's of arguments to runDMCStep
   *
   *  There should be a division between const input to runVMCStep
   *  And step to step state
   */
  struct StateForThread
  {
    const QMCDriverInput& qmcdrv_input;
    const DriftModifierBase& drift_modifier;
    const MCPopulation& population;
    SFNBranch& branch_engine;
    IndexType recalculate_properties_period;
    const size_t steps_per_block;
    IndexType step            = -1;
    IndexType global_step     = -1;
    bool is_recomputing_block = false;
    /// if true, calculating walker one-by-one within a crowd
    const bool serializing_crowd_walkers;

    StateForThread(const QMCDriverInput& qmci,
                   DriftModifierBase& drift_mod,
                   SFNBranch& branch_eng,
                   MCPopulation& pop,
                   const size_t steps_per_block,
                   const bool serializing_crowd_walkers)
        : qmcdrv_input(qmci),
          drift_modifier(drift_mod),
          population(pop),
          branch_engine(branch_eng),
          steps_per_block(steps_per_block),
          serializing_crowd_walkers(serializing_crowd_walkers)
    {}
  };

  class DMCTimers
  {
  public:
    NewTimer& tmove_timer;
    NewTimer& step_begin_recompute_timer;
    DMCTimers(const std::string& prefix)
        : tmove_timer(createGlobalTimer(prefix + "Tmove", timer_level_medium)),
          step_begin_recompute_timer(createGlobalTimer(prefix + "Step_begin_recompute", timer_level_medium))
    {}
  };

  /// Constructor.
  DMCBatched(const ProjectData& project_data,
             QMCDriverInput&& qmcdriver_input,
             UPtr<EstimatorManagerNew>&& estimator_manager,
             DMCDriverInput&& input,
             WalkerConfigurations& wc,
             MCPopulation&& pop,
             const RefVector<RandomBase<FullPrecRealType>>& rng_refs,
             Communicate* comm);

  /// Copy Constructor (disabled)
  DMCBatched(const DMCBatched&) = delete;
  /// Copy operator (disabled).
  DMCBatched& operator=(const DMCBatched&) = delete;

  DMCBatched(DMCBatched&&) = default;

  ~DMCBatched() override;

  /** DMCBatched driver will eventually ignore cur
   *
   *  This is the shared entry point
   *  from QMCMain so cannot be updated yet
   *
   *  Contains logic that sets walkers_per_rank_
   *  TargetWalkers trump walkers, if it is not set
   *  walkers which is by default per rank for the batched drivers
   *  from this or the previous section wins.
   *
   *  walkers is still badly named.
   */
  void process(xmlNodePtr cur) override;

  bool run() override;

  QMCRunType getRunType() override { return QMCRunType::DMC_BATCH; }

private:
  /// forward declaration. DMC specialized ContextForSteps
  class DMCContextForSteps;

  const DMCDriverInput dmcdriver_input_;
  /// Per crowd, driver-specific move contexts
  UPtrVector<DMCContextForSteps> step_contexts_;
  /// obtain reference vector of step contexts
  RefVector<ContextForSteps> getContextForStepsRefs() const;
  /** I think its better if these have there own type and variable name
   */
  DMCTimers dmc_timers_;
  /// Interval between branching
  IndexType branch_interval_;
  ///branch engine
  std::unique_ptr<SFNBranch> branch_engine_;
  ///walker controller for load-balance
  std::unique_ptr<WalkerControl> walker_controller_;

  // create Rngs and StepContests
  void createStepContexts(int num_crowds);

  struct RunContext
  {
    WalkerLogManager wlog_manager;
    LoopTimer<> dmc_loop;
    RunTimeControl<> runtimeControl;
    StateForThread dmc_state;
  };

  RunContext startRun();

  // This is the task body executed at crowd scope
  // it does not have access to object members by design
  static void runDMCStep(int crowd_id,
                         const StateForThread& sft,
                         DriverTimers& timers,
                         DMCTimers& dmc_timers,
                         UPtrVector<DMCContextForSteps>& move_context,
                         UPtrVector<Crowd>& crowds);

  template<CoordsType CT>
  static void advanceWalkers(const StateForThread& sft,
                             Crowd& crowd,
                             DriverTimers& timers,
                             DMCTimers& dmc_timers,
                             DMCContextForSteps& move_context,
                             bool recompute,
                             bool accumulate_this_step);

  /** @ingroup testing helper functions
   * @brief integration testing requires an ability to do just some of
   * what the application functions demand
   * @{
   */
  StateForThread mockRunStart();
  auto getCrowds() -> decltype(crowds_)& { return crowds_; }
  /// @}
  friend class qmcplusplus::testing::DMCBatchedTest;
  friend class qmcplusplus::testing::DMCBatchedTestAccessor;
};

} // namespace qmcplusplus

#endif
