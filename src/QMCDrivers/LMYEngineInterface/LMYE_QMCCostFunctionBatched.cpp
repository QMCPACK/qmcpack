//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu, University of California, Berkeley
//
// File created by: Leon Otis, leon_otis@berkeley.edu, University of California, Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/WFOpt/QMCCostFunctionBatched.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/CommOperators.h"

#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"

namespace qmcplusplus
{
size_t QMCCostFunctionBatched::total_samples() { return samples_.getGlobalNumSamples(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the cost function using the LMYEngine for interfacing with batched driver
///
///////////////////////////////////////////////////////////////////////////////////////////////////
QMCCostFunctionBatched::Return_rt QMCCostFunctionBatched::LMYEngineCost_detail(
    cqmc::engine::LMYEngine<Return_t>* EngineObj)
{
  // get total number of samples
  const size_t m = this->total_samples();
  // reset Engine object
  EngineObj->reset();

  // turn off wavefunction update mode
  EngineObj->turn_off_update();

#pragma omp parallel
  {
    int ip = omp_get_thread_num();
    // for each thread, loop over samples
    const int nw = rank_local_num_samples_;
    for (int iw = 0; iw < nw; iw++)
    {
      // get a pointer to the record for this sample
      const Return_rt* restrict saved = RecordsOnNode_[iw];

      // take sample
      EngineObj->take_sample(saved[ENERGY_NEW], 1.0, saved[REWEIGHT] / SumValue[SUM_WGT]);
    }
  }
  //}
  // finish taking sample
  EngineObj->sample_finish();

  // compute energy and target relevant quantities
  EngineObj->energy_target_compute();

  // prepare variables to hold the output of the engine call
  double energy_avg  = EngineObj->energy_mean();
  double energy_sdev = EngineObj->energy_sdev();
  double energy_serr = EngineObj->energy_statistical_err();
  double target_avg  = EngineObj->target_value();
  double target_serr = EngineObj->target_statistical_err();


  // return the cost function value (target function if we are targeting excited states and energy if we are doing groud state calculations)
  double cost_value = (targetExcited ? target_avg : energy_avg);
  return cost_value;
}

} // namespace qmcplusplus
