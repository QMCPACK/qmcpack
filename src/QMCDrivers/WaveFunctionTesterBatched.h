//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WAVEFUNCTIONTESTERBATCHED_H
#define QMCPLUSPLUS_WAVEFUNCTIONTESTERBATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/WaveFunctionTesterBatchedInput.h"

namespace qmcplusplus
{
class WaveFunctionTesterBatched : public QMCDriverNew
{
public:
  /** Construct a crowd-native wavefunction tester with its parsed diagnostic settings. */
  WaveFunctionTesterBatched(const ProjectData& project_data,
                            QMCDriverInput&& qmcdriver_input,
                            UPtr<EstimatorManagerNew>&& estimator_manager,
                            WaveFunctionTesterBatchedInput&& input,
                            WalkerConfigurations& wc,
                            MCPopulation&& population,
                            const RefVector<RandomBase<FullPrecRealType>>& rng_refs,
                            Communicate* comm);

  /** Create the diagnostic walker population, crowds, and per-crowd random contexts. */
  void process(xmlNodePtr node) override;
  /** Evaluate finite-difference and ratio diagnostics and write rank-local wftest output. */
  bool run() override;

private:
  class TestContext : public ContextForSteps
  {
  public:
    /** Bind the crowd diagnostic context to its assigned random-number generator. */
    explicit TestContext(RandomBase<FullPrecRealType>& random_gen) : ContextForSteps(random_gen) {}
  };

  WaveFunctionTesterBatchedInput input_;
  UPtrVector<ContextForSteps> test_contexts_;

  /** Identify this driver as the batched wavefunction-test implementation. */
  QMCRunType getRunType() override { return QMCRunType::WF_TEST_BATCH; }
};
} // namespace qmcplusplus

#endif