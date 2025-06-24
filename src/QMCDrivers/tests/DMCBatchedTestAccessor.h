
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DMCBATCHEDTESTACCESSOR_H
#define QMCPLUSPLUS_DMCBATCHEDTESTACCESSOR_H

#include "DMC/DMCBatched.h"

namespace qmcplusplus::testing
{
class DMCBatchedTestAccessor
{
public:
  using DMCContextForSteps = DMCBatched::DMCContextForSteps;
  using DMCTimers          = DMCBatched::DMCTimers;
  using StateForThread     = DMCBatched::StateForThread;
  using DriverTimers       = DMCBatched::DriverTimers;
  static auto mockRunStart(DMCBatched& dmcb) -> decltype(dmcb.mockRunStart()) { return dmcb.mockRunStart(); }
  static auto getContextForSteps(DMCBatched& dmcb) -> decltype(dmcb.step_contexts_)& { return dmcb.step_contexts_; }
  static auto getCrowds(DMCBatched& dmcb) -> decltype(dmcb.getCrowds()) { return dmcb.getCrowds(); }
  static std::size_t getStepsPerBlock(DMCBatched& dmcb) { return dmcb.steps_per_block_; }
  static auto startRun(DMCBatched& dmcb) -> decltype(dmcb.startRun()) { return dmcb.startRun(); }
  static auto getTimers(DMCBatched& dmcb) -> decltype(dmcb.timers_)& { return dmcb.timers_; }
  static auto getDMCTimers(DMCBatched& dmcb) -> decltype(dmcb.dmc_timers_)& { return dmcb.dmc_timers_; }
  static void runDMCStep(DMCBatched& dmcb,
                         int crowd_id,
                         const StateForThread& sft,
                         DriverTimers& timers,
                         DMCTimers& dmc_timers,
                         UPtrVector<DMCContextForSteps>& step_context,
                         UPtrVector<Crowd>& crowds)
  {
    DMCBatched::runDMCStep(crowd_id, sft, timers, dmc_timers, step_context, crowds);
  }
  static void endBlock(DMCBatched& dmcb, int block)
  {
    dmcb.endBlock();
    dmcb.recordBlock(block);
  }
  static void endRun(DMCBatched& dmcb, int block)
  {
    dmcb.estimator_manager_->stopDriverRun();
    dmcb.finalize(block);
  }
};

} // namespace qmcplusplus::testing

#endif
