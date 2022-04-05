//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COSTFUNCTIONCROWDDATA_H
#define QMCPLUSPLUS_COSTFUNCTIONCROWDDATA_H

#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "DriverWalkerTypes.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 */


/// Class to hold temporary data and object copies for crowd-local evaluation
class CostFunctionCrowdData
{
public:
  using Return_rt = qmcplusplus::QMCTraits::RealType;

  /// Create the arrays of crowd_size and create object copies
  CostFunctionCrowdData(int crowd_size,
                        ParticleSet& P,
                        TrialWaveFunction& Psi,
                        QMCHamiltonian& H,
                        std::vector<std::string>& H_KE_node_names,
                        RandomGenerator& Rng);

  /// Set the log_psi_* arrays to zero
  void zero_log_psi();

  RefVector<ParticleSet> get_p_list(int len);
  RefVector<TrialWaveFunction> get_wf_list(int len);
  RefVector<QMCHamiltonian> get_h_list(int len);
  RefVector<QMCHamiltonian> get_h0_list(int len);

  std::vector<Return_rt>& get_log_psi_fixed() { return log_psi_fixed_; }
  std::vector<Return_rt>& get_log_psi_opt() { return log_psi_opt_; }

  UPtrVector<RandomGenerator>& get_rng_ptr_list() { return rng_ptr_list_; }
  RandomGenerator& get_rng_save() { return *rng_save_ptr_; }

  UPtrVector<TrialWaveFunction>& get_wf_ptr_list() { return wf_ptr_list_; }

  Return_rt& get_e0() { return e0_; }
  Return_rt& get_e2() { return e2_; }

  Return_rt& get_wgt() { return wgt_; }
  Return_rt& get_wgt2() { return wgt2_; }

  DriverWalkerResourceCollection& getSharedResource() { return driverwalker_resource_collection_; }

private:
  // Temporary vectors for the call to flex_evaluateDeltaLogSetup
  std::vector<Return_rt> log_psi_fixed_;
  std::vector<Return_rt> log_psi_opt_;

  // List of objects for use in flex_* calls
  UPtrVector<TrialWaveFunction> wf_ptr_list_;
  UPtrVector<ParticleSet> p_ptr_list_;
  UPtrVector<QMCHamiltonian> h_ptr_list_;
  UPtrVector<QMCHamiltonian> h0_ptr_list_;
  UPtrVector<RandomGenerator> rng_ptr_list_;

  // proivides multi walker resource
  DriverWalkerResourceCollection driverwalker_resource_collection_;

  // Saved RNG state to reset to before correlated sampling
  std::unique_ptr<RandomGenerator> rng_save_ptr_;

  // Crowd-local accumulator variables
  Return_rt e0_;
  Return_rt e2_;

  Return_rt wgt_;
  Return_rt wgt2_;
};


} // namespace qmcplusplus
#endif
