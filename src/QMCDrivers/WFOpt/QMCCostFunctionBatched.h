//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COSTFUNCTION_H
#define QMCPLUSPLUS_COSTFUNCTION_H

#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC running on a single thread.
 */


/// Class to hold temporary data and object copies for crowd-local evaluation
class OptimizerEvalData
{
public:
  typedef qmcplusplus::QMCTraits::RealType Return_rt;

  /// Create the arrays of crowd_size and create object copies
  OptimizerEvalData(int crowd_size,
                    ParticleSet& P,
                    TrialWaveFunction& Psi,
                    QMCHamiltonian& H,
                    QMCHamiltonian& H_KE_Node,
                    RandomGenerator_t& Rng);

  /// Creates the reference vectors.  On the last iteration these may be smaller than the full crowd size
  void set_current_crowd_size(int curr_crowd_size);

  /// Set the log_psi_* arrays to zero
  void zero_log_psi();

  RefVector<ParticleSet>& p_list() { return p_list_; }
  RefVector<TrialWaveFunction>& wf_list() { return wf_list_; }
  RefVector<QMCHamiltonian>& h_list() { return h_list_; }
  RefVector<QMCHamiltonian>& h0_list() { return h0_list_; }

  std::vector<Return_rt>& log_psi_fixed() { return log_psi_fixed_; }
  std::vector<Return_rt>& log_psi_opt() { return log_psi_opt_; }

  UPtrVector<RandomGenerator_t>& get_rng_ptr_list() { return rng_ptr_list_; }
  RandomGenerator_t& get_rng_save() { return *rng_save_ptr_; }

  UPtrVector<TrialWaveFunction>& wf_ptr_list() { return wf_ptr_list_; }

  Return_rt& e0() { return e0_; }
  Return_rt& e2() { return e2_; }

  Return_rt& wgt() { return wgt_; }
  Return_rt& wgt2() { return wgt2_; }


private:
  int crowd_size_;
  int crowd_size_curr_;

  std::vector<Return_rt> log_psi_fixed_;
  std::vector<Return_rt> log_psi_opt_;

  UPtrVector<TrialWaveFunction> wf_ptr_list_;
  UPtrVector<ParticleSet> p_ptr_list_;
  UPtrVector<QMCHamiltonian> h_ptr_list_;
  UPtrVector<QMCHamiltonian> h0_ptr_list_;
  UPtrVector<RandomGenerator_t> rng_ptr_list_;
  std::unique_ptr<RandomGenerator_t> rng_save_ptr_;

  RefVector<TrialWaveFunction> wf_list_;
  RefVector<ParticleSet> p_list_;
  RefVector<QMCHamiltonian> h_list_;
  RefVector<QMCHamiltonian> h0_list_;

  // Crowd-local accumulator variables
  Return_rt e0_;
  Return_rt e2_;

  Return_rt wgt_;
  Return_rt wgt2_;
};


class QMCCostFunctionBatched : public QMCCostFunctionBase, public QMCTraits
{
public:
  ///Constructor.
  QMCCostFunctionBatched(MCWalkerConfiguration& w,
                         TrialWaveFunction& psi,
                         QMCHamiltonian& h,
                         SampleStack& samples,
                         int opt_num_crowds,
                         int crowd_size,
                         Communicate* comm);

  ///Destructor
  ~QMCCostFunctionBatched();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
#ifdef HAVE_LMY_ENGINE
  void engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                  DescentEngine& descentEngineObj,
                                  const std::string& MinMethod);
#endif


  void resetPsi(bool final_reset = false);
  void GradCost(std::vector<Return_rt>& PGradient, const std::vector<Return_rt>& PM, Return_rt FiniteDiff = 0);
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right);

protected:
  std::unique_ptr<QMCHamiltonian> H_KE_Node;
  std::unique_ptr<QMCHamiltonian> extractFixedHamiltonianComponents();

  Matrix<Return_rt> RecordsOnNode_;

  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  Matrix<Return_rt> DerivRecords_;
  Matrix<Return_rt> HDerivRecords_;

  Return_rt correlatedSampling(bool needGrad = true);

  SampleStack& samples_;

  int opt_batch_size_;
  int opt_num_crowds_;

  std::vector<std::unique_ptr<OptimizerEvalData>> opt_eval_;

  NewTimer& check_config_timer_;
  NewTimer& corr_sampling_timer_;
  NewTimer& fill_timer_;


#ifdef HAVE_LMY_ENGINE
  int total_samples();
#endif
};
} // namespace qmcplusplus
#endif
