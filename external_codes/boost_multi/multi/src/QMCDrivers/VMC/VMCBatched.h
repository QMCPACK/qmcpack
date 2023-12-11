//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_VMCBATCHED_H
#define QMCPLUSPLUS_VMCBATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/ContextForSteps.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"

#include "Utilities/Timer.h"

namespace qmcplusplus
{
namespace testing
{
class VMCBatchedTest;
}

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCBatched : public QMCDriverNew
{
public:
  using Base              = QMCDriverNew;
  using FullPrecRealType  = QMCTraits::FullPrecRealType;
  using PosType           = QMCTraits::PosType;
  using ParticlePositions = PtclOnLatticeTraits::ParticlePos;

  /** To avoid 10's of arguments to runVMCStep
   *
   *  There should be a division between const input to runVMCStep
   *  And step to step state
   */
  struct StateForThread
  {
    const QMCDriverInput& qmcdrv_input;
    const VMCDriverInput& vmcdrv_input;
    const DriftModifierBase& drift_modifier;
    const MCPopulation& population;
    IndexType recalculate_properties_period;
    IndexType step            = -1;
    bool is_recomputing_block = false;

    StateForThread(const QMCDriverInput& qmci,
                   const VMCDriverInput& vmci,
                   DriftModifierBase& drift_mod,
                   MCPopulation& pop)
        : qmcdrv_input(qmci), vmcdrv_input(vmci), drift_modifier(drift_mod), population(pop)
    {}
  };

public:
  /// Constructor.
  VMCBatched(const ProjectData& project_data,
             QMCDriverInput&& qmcdriver_input,
             const std::optional<EstimatorManagerInput>& global_emi,
             VMCDriverInput&& input,
             WalkerConfigurations& wc,
             MCPopulation&& pop,
             SampleStack& samples_,
             Communicate* comm);

  void process(xmlNodePtr node) override;

  bool run() override;

  /** Refactor of VMCUpdatePbyP in crowd context
   *
   *  MCWalkerConfiguration layer removed.
   *  Obfuscation of state changes via buffer and MCWalkerconfiguration require this be tested well
   */
  template<CoordsType CT>
  static void advanceWalkers(const StateForThread& sft,
                             Crowd& crowd,
                             DriverTimers& timers,
                             ContextForSteps& move_context,
                             bool recompute,
                             bool accumulate_this_step);

  // This is the task body executed at crowd scope
  // it does not have access to object member variables by design
  static void runVMCStep(int crowd_id,
                         const StateForThread& sft,
                         DriverTimers& timers,
                         UPtrVector<ContextForSteps>& context_for_steps,
                         UPtrVector<Crowd>& crowds);

  /** transitional interface on the way to better walker count adjustment handling.
   *  returns a closure taking walkers per rank and accomplishing what calc_default_local_walkers does.
   */
  auto getCDLW();

  /** Enable collecting samples during the VMC run
   *
   *  strong assumption that VMCBatched driver has passed through process phase of
   *  initialization.
   *  A side effect of VMCBatched::process is that MCPopulation has created local walkers.
   */
  void enable_sample_collection();

private:
  int prevSteps;
  int prevStepsBetweenSamples;
  VMCDriverInput vmcdriver_input_;
  QMCRunType getRunType() override { return QMCRunType::VMC_BATCH; }
  /// Copy constructor
  VMCBatched(const VMCBatched&) = delete;
  /// Copy operator (disabled).
  VMCBatched& operator=(const VMCBatched&) = delete;


  /// Storage for samples (later used in optimizer)
  SampleStack& samples_;
  /// Sample collection flag
  bool collect_samples_;
  /** function to calculate samples per MPI rank
   */
  static int compute_samples_per_rank(const QMCDriverInput& qmcdriver_input, const IndexType local_walkers);

  friend class qmcplusplus::testing::VMCBatchedTest;
};

extern std::ostream& operator<<(std::ostream& o_stream, const VMCBatched& vmc_batched);
} // namespace qmcplusplus

#endif
