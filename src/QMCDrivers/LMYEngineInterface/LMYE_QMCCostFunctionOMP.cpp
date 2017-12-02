//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luning Zhao, zhaoln@berkeley.edu, University of California, Berkeley
//                    Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//
// File created by: Luning Zhao, zhaoln@berkeley.edu, University of California, Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/QMCCostFunctionOMP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"

//#include "Eigen/Dense"
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"

namespace qmcplusplus {

int QMCCostFunctionOMP::total_samples() {

  // for the unfamiliar, the [] starts a lambda function
  return std::accumulate(wClones.begin(), wClones.begin()+NumThreads, 0, [] (int x, const MCWalkerConfiguration * p) { return x + p->numSamples(); });

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Computes the cost function using the LMYEngine
///
///////////////////////////////////////////////////////////////////////////////////////////////////
QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::LMYEngineCost_detail(cqmc::engine::LMYEngine * EngineObj) {

  // get total number of samples
  const int m = this->total_samples();

  // reset Engine object 
  EngineObj->reset();

  // turn off wavefunction update mode 
  EngineObj->turn_off_update();

  //for (int ip = 0, j = 0; ip < NumThreads; ip++) {
# pragma omp parallel
  {
    int ip = omp_get_thread_num();     
    // for each thread, loop over samples
    const int nw = wClones[ip]->numSamples();
    for (int iw = 0; iw < nw; iw++) {

      // get a pointer to the record for this sample
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];

      // record this sample's weight (normalized by division by the sum of weights I think...)
      //wgt_vec.at(j) = saved[REWEIGHT] / SumValue[SUM_WGT];

      // record the sample's local energy
      //lce_vec.at(j) = saved[ENERGY_NEW];

      // take sample 
      EngineObj->take_sample(saved[ENERGY_NEW], 1.0, saved[REWEIGHT] / SumValue[SUM_WGT]);

    }
  }
  //}
  // finish taking sample 
  EngineObj->sample_finish();

  // compute energy and target relevant quantities 
  EngineObj->energy_target_compute();

  // prepare variables to hold the ouput of the engine call
  double energy_avg  = EngineObj->energy_mean(); 
  double energy_sdev = EngineObj->energy_sdev(); 
  double energy_serr = EngineObj->energy_statistical_err();  
  double target_avg  = EngineObj->target_value();
  double target_serr = EngineObj->target_statistical_err();

  // prepare a stream to hold engine printout
  //std::stringstream engine_out;

  // call the LMYEngine to compute the cost function
  //cqmc::engine::call_engine(false, // exact_sampling
  //                          !targetExcited, // ground_state
  //                          false, // variance_correct
  //                          true, // print
  //                          omega_shift, // hd_lm_shift
  //                          0.0, // var_weight
  //                          lce_vec, vgs_vec, wgt_vec,
  //                          energy_avg, energy_sdev, energy_serr,
  //                          target_avg, target_serr, engine_out);

  // print the engine output
  //app_log() << engine_out.str();

  // return the cost function value (target function if we are targeting excited states and energy if we are doing groud state calculations)
  double cost_value = (targetExcited ? target_avg : energy_avg);
  return cost_value;

}

}
