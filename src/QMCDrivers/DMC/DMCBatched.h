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
#include "QMCDrivers/ContextForSteps.h"
#include "Particle/MCCoords.hpp"

namespace qmcplusplus
{
class DriverModifierBase;
class WalkerControl;
class SFNBranch;

namespace testing
{
class DMCBatchedTest;
}
/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a DMC using particle-by-particle threaded and batched moves.
 */
class DMCBatched : public QMCDriverNew
{
public:
  using Base              = QMCDriverNew;
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
    const DMCDriverInput& dmcdrv_input;
    const DriftModifierBase& drift_modifier;
    const MCPopulation& population;
    SFNBranch& branch_engine;
    IndexType recalculate_properties_period;
    IndexType step            = -1;
    bool is_recomputing_block = false;
    StateForThread(const QMCDriverInput& qmci,
                   const DMCDriverInput& dmci,
                   DriftModifierBase& drift_mod,
                   SFNBranch& branch_eng,
                   MCPopulation& pop)
        : qmcdrv_input(qmci), dmcdrv_input(dmci), drift_modifier(drift_mod), population(pop), branch_engine(branch_eng)
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
             const std::optional<EstimatorManagerInput>& global_emi,
             DMCDriverInput&& input,
             WalkerConfigurations& wc,
             MCPopulation&& pop,
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

  // This is the task body executed at crowd scope
  // it does not have access to object members by design
  static void runDMCStep(int crowd_id,
                         const StateForThread& sft,
                         DriverTimers& timers,
                         DMCTimers& dmc_timers,
                         UPtrVector<ContextForSteps>& move_context,
                         UPtrVector<Crowd>& crowds);


  QMCRunType getRunType() override { return QMCRunType::DMC_BATCH; }

  void setNonLocalMoveHandler(QMCHamiltonian& golden_hamiltonian);

private:
  const DMCDriverInput dmcdriver_input_;

  /** I think its better if these have there own type and variable name
   */
  DMCTimers dmc_timers_;
  /// Interval between branching
  IndexType branch_interval_;
  ///branch engine
  std::unique_ptr<SFNBranch> branch_engine_;
  ///walker controller for load-balance
  std::unique_ptr<WalkerControl> walker_controller_;

  template<CoordsType CT>
  static void advanceWalkers(const StateForThread& sft,
                             Crowd& crowd,
                             DriverTimers& timers,
                             DMCTimers& dmc_timers,
                             ContextForSteps& move_context,
                             bool recompute,
                             bool accumulate_this_step);

  friend class qmcplusplus::testing::DMCBatchedTest;
};

} // namespace qmcplusplus

#endif
