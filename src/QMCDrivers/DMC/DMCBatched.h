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
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos_t;
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
    IndexType step = -1;
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
    DMCTimers(const std::string& prefix) : tmove_timer(*timer_manager.createTimer(prefix + "Tmove", timer_level_medium))
    {}
  };

  /// Constructor.
  DMCBatched(const ProjectData& project_data,
             QMCDriverInput&& qmcdriver_input,
             DMCDriverInput&& input,
             MCPopulation&& pop,
             TrialWaveFunction& psi,
             QMCHamiltonian& h,
             Communicate* comm);

  /// Copy Constructor (disabled)
  DMCBatched(const DMCBatched&) = delete;
  /// Copy operator (disabled).
  DMCBatched& operator=(const DMCBatched&) = delete;

  DMCBatched(DMCBatched&&) = default;

  ~DMCBatched();

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
  void process(xmlNodePtr cur);

  bool run();

  // This is the task body executed at crowd scope
  // it does not have access to object members by design
  static void runDMCStep(int crowd_id,
                         const StateForThread& sft,
                         DriverTimers& timers,
                         DMCTimers& dmc_timers,
                         UPtrVector<ContextForSteps>& move_context,
                         UPtrVector<Crowd>& crowds);


  QMCRunType getRunType() { return QMCRunType::DMC_BATCH; }

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

  static void advanceWalkers(const StateForThread& sft,
                             Crowd& crowd,
                             DriverTimers& timers,
                             DMCTimers& dmc_timers,
                             ContextForSteps& move_context,
                             bool recompute);

  /** Allows us to build complete reference set for walkers.
   */
  struct DMCPerWalkerRefs
  {
    RefVector<MCPWalker> walkers;
    RefVector<TrialWaveFunction> walker_twfs;
    RefVector<QMCHamiltonian> walker_hamiltonians;
    RefVector<ParticleSet> walker_elecs;
    RefVector<FullPrecRealType> old_energies;
    RefVector<FullPrecRealType> new_energies;
    RefVector<RealType> rr_proposed;
    RefVector<RealType> rr_accepted;
    RefVector<RealType> gf_accs;
    DMCPerWalkerRefs(int nwalkers)
    {
      walkers.reserve(nwalkers);
      walker_twfs.reserve(nwalkers);
      walker_hamiltonians.reserve(nwalkers);
      walker_elecs.reserve(nwalkers);
      old_energies.reserve(nwalkers);
      new_energies.reserve(nwalkers);
      rr_proposed.reserve(nwalkers);
      rr_accepted.reserve(nwalkers);
      gf_accs.reserve(nwalkers);
    }
  };

  /** To allow breaking up implementation without creating giant function signature.
   */
  struct DMCPerWalkerRefRefs
  {
    const RefVector<MCPWalker>& walkers;
    const RefVector<TrialWaveFunction>& walker_twfs;
    const RefVector<QMCHamiltonian>& walker_hamiltonians;
    const RefVector<ParticleSet>& walker_elecs;
    std::vector<FullPrecRealType>& old_energies;
    std::vector<FullPrecRealType>& new_energies;
    std::vector<RealType>& rr_proposed;
    std::vector<RealType>& rr_accepted;
    std::vector<RealType>& gf_accs;
  };


  /** for the return of DMCPerWalkerRefs split into moved and stalled
   *
   *  until C++17 we need a structure to return the split moved and stalled refs
   *  I find this more readable than the tuple tie pattern
   */
  struct MovedStalled
  {
    MovedStalled(int num_walkers, int num_moved) : moved(num_moved), stalled(num_walkers - num_moved) {}
    DMCPerWalkerRefs moved;
    DMCPerWalkerRefs stalled;
  };

  static MovedStalled buildMovedStalled(const std::vector<int>& did_walker_move, const DMCPerWalkerRefRefs& refs);

  static void handleMovedWalkers(DMCPerWalkerRefs& moved,
                                 const StateForThread& sft,
                                 DriverTimers& timers,
                                 bool recompute);
  static void handleStalledWalkers(DMCPerWalkerRefs& stalled, const StateForThread& sft, bool recompute);
  // struct DMCTimers
  // {
  //   NewTimer& dmc_movePbyP;
  //   DriverTimers(const std::string& prefix)
  //       : dmc_movePbyP(*timer_manager.createTimer(prefix + "DMC_movePbyP", timer_level_medium)),
  //   {}
  // };

  // DMCTimers dmc_timers_;

  friend class qmcplusplus::testing::DMCBatchedTest;
};

} // namespace qmcplusplus

#endif
