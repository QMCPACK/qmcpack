//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu University, University
// of California Berkeley
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Leon Otis, leon_otis@berkeley.edu University, University of
// California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DESCENT_ENGINE_HEADER
#define QMCPLUSPLUS_DESCENT_ENGINE_HEADER

#include "Configuration.h"
#include "Message/Communicate.h"
#include "VariableSet.h"
#include <libxml/tree.h>
#include <vector>

namespace qmcplusplus
{
class DescentEngine
{
  using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
  using ValueType         = qmcplusplus::QMCTraits::ValueType;
  using RealType          = qmcplusplus::QMCTraits::RealType;
  using FullPrecRealType  = qmcplusplus::QMCTraits::FullPrecRealType;

private:
  ///Vector for local energy parameter derivatives
  std::vector<FullPrecValueType> avg_le_der_samp_;
  ///Vector for local energy parameter derivatives on one thread
  std::vector<std::vector<FullPrecValueType>> replica_le_der_samp_;

  ///Vector for WF parameter derivatives
  std::vector<FullPrecValueType> avg_der_rat_samp_;
  ///Vector for WF parameter derivatives on one thread
  std::vector<std::vector<FullPrecValueType>> replica_der_rat_samp_;

  ///Vector for target function numerator parameter derivatives
  std::vector<FullPrecValueType> avg_numer_der_samp_;
  ///Vector for target function numerator parameter derivatives on one thread
  std::vector<std::vector<FullPrecValueType>> replica_numer_der_samp_;

  ///Vector for target function denominator parameter derivatives
  std::vector<FullPrecValueType> avg_denom_der_samp_;
  ///Vector for target function denominator parameter derivatives on one thread
  std::vector<std::vector<FullPrecValueType>> replica_denom_der_samp_;

  ///Total sum of weights
  ValueType w_sum_;

  ///Average energy on a descent step
  ValueType e_avg_;
  ///Variance of the energy
  ValueType e_var_;
  ///Standard deviation of the energy
  ValueType e_sd_;
  ///Standard error of the energy
  ValueType e_err_;

  ///Average target function numerator on a descent step
  ValueType numer_avg_;
  ///Variance of the target function numerator
  ValueType numer_var_;
  ///Standard error of the target function numerator
  ValueType numer_err_;

  ///Average target function denominator on a descent step
  ValueType denom_avg_;
  ///Variance of the target function denominator
  ValueType denom_var_;
  ///Standard error of the target function denominator
  ValueType denom_err_;

  ///Average target function value on a descent step
  ValueType target_avg_;
  ///Variance of the target function
  ValueType target_var_;
  ///Standard error of the target function
  ValueType target_err_;

  /// history of sampled |value/guiding|^2 ratios for one iteration
  std::vector<ValueType> vg_history_;
  std::vector<std::vector<ValueType>> replica_vg_history_;

  /// history of sampled |value/guiding|^2 ratios  during the descent
  /// finalization phase
  std::vector<ValueType> final_vg_history_;
  std::vector<std::vector<ValueType>> replica_final_vg_history_;


  /// history of sampled configuration weights for one iteration
  std::vector<ValueType> w_history_;
  std::vector<std::vector<ValueType>> replica_w_history_;
  /// history of sampled configuration weights during descent finalization phase
  std::vector<ValueType> final_w_history_;
  std::vector<std::vector<ValueType>> replica_final_w_history_;

  /// a history of sampled local energies times the |value/guiding|^2 raitos for
  /// one iteration
  std::vector<ValueType> lev_history_;
  std::vector<std::vector<ValueType>> replica_lev_history_;

  /// history of sampled local energies times the |value/guiding|^2 raitos during
  /// the descent finalization phase
  std::vector<ValueType> final_lev_history_;
  std::vector<std::vector<ValueType>> replica_final_lev_history_;

  /// a vector to store the averages of the energy during the descent
  /// finalization phase
  std::vector<ValueType> final_le_avg_history_;

  /// a vector to store the variances of the energy during the descent
  /// finalization phase
  std::vector<ValueType> final_var_avg_history_;

  /// a history of target function numerator times the |value/guiding|^2 ratios
  /// for one iteration
  std::vector<ValueType> tnv_history_;
  std::vector<std::vector<ValueType>> replica_tnv_history_;

  /// a history of target function numerator times the |value/guiding|^2 ratios
  /// during the descent finalization phase
  std::vector<ValueType> final_tnv_history_;
  std::vector<std::vector<ValueType>> replica_final_tnv_history_;

  /// a history of target function denominator times the |value/guiding|^2 ratios
  /// for one iteration
  std::vector<ValueType> tdv_history_;
  std::vector<std::vector<ValueType>> replica_tdv_history_;
  /// a history of target function denomerator times the |value/guiding|^2 ratios
  /// during the descent finalization phase
  std::vector<ValueType> final_tdv_history_;
  std::vector<std::vector<ValueType>> replica_final_tdv_history_;

  /// a vector to store the averages of the target function during the descent
  /// finalization phase
  std::vector<ValueType> final_tar_avg_history_;
  /// a vector to store the variances of the target function during the descent
  /// finalization phase
  std::vector<ValueType> final_tar_var_history_;

  /// Vector that stores the final averaged derivatives of the cost function
  std::vector<ValueType> lderivs_;

  /// Communicator handles MPI reduction
  Communicate* my_comm_;

  /// Whether to target excited state
  bool engine_target_excited_;

  /// Number of optimizable parameters
  int num_params_;

  /// Vector for storing parameter values from previous optimization step
  std::vector<ValueType> params_copy_;

  /// Vector for storing parameter values for current optimization step
  std::vector<ValueType> current_params_;

  /// Vector for storing Lagrangian derivatives from previous optimization steps
  std::vector<std::vector<ValueType>> deriv_records_;

  /// Vector for storing step size denominator values from previous optimization
  /// step
  std::vector<ValueType> denom_records_;

  /// Vector for storing step size numerator values from previous optimization
  /// step
  std::vector<ValueType> numer_records_;

  /// Parameter for accelerated descent recursion relation
  ValueType lambda_ = 0.0;
  /// Vector for storing step sizes from previous optimization step.
  std::vector<ValueType> taus_;
  /// Vector for storing running average of squares of the derivatives
  std::vector<ValueType> derivs_squared_;

  /// Integer for keeping track of only number of descent steps taken
  int descent_num_;

  /// What variety of gradient descent will be used
  std::string flavor_;

  /// Step sizes for different types of parameters
  ValueType tjf_2body_eta_;
  ValueType tjf_1body_eta_;
  ValueType f_eta_;
  ValueType gauss_eta_;
  ValueType ci_eta_;
  ValueType orb_eta_;

  /// Whether to gradually ramp up step sizes in descent
  bool ramp_eta_;

  /// Number of steps over which to ramp up step size
  int ramp_num_;

  /// Number of parameter difference vectors stored when descent is used in a
  /// hybrid optimization
  int store_num_;

  /// Counter of how many vectors have been stored so far
  int store_count_;

  /// Vectors of parameter names and types, used in the assignment of step sizes
  std::vector<std::string> engine_param_names_;
  std::vector<int> engine_param_types_;

  /// Vector for storing parameter values for calculating differences to be given
  /// to hybrid method
  std::vector<ValueType> params_for_diff_;

  /// Vector for storing the input vectors to the BLM steps of hybrid method
  std::vector<std::vector<ValueType>> hybrid_blm_input_;

  /// Value of omega in excited state functional
  ValueType omega_;

  /// Iteration to start collecting samples for final average and error blocking
  /// analysis
  int collection_step_;

  /// Iteration to start computing averages and errors from the stored values
  /// during the finalization phase
  int compute_step_;

  /// Whether to start collecting samples for the histories in the finalization
  /// phase
  bool collect_count_ = false;

  /// Counter for the number of descent steps taken in the finalization phase
  int final_descent_num_ = 0;

  /// Whether to print out derivative terms for each parameter
  std::string print_deriv_;

public:
  /// Constructor for engine
  DescentEngine(Communicate* comm, const xmlNodePtr cur);

  /// process xml node
  bool processXML(const xmlNodePtr cur);

  /// Prepare for taking samples
  void prepareStorage(const int num_replicas, const int num_optimizables);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that Take Sample Data from the Host Code
  ///
  /// \param[in]  der_rat_samp   <n|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var
  /// )
  /// \param[in]  le_der_samp    <n|H|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ...
  /// N_var )
  /// \param[in]  ls_der_samp    <|S^2|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ...
  /// N_var )
  /// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
  /// \param[in]  weight_samp    weight for this sample
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void takeSample(const int replica_id,
                  const std::vector<FullPrecValueType>& der_rat_samp,
                  const std::vector<FullPrecValueType>& le_der_samp,
                  const std::vector<FullPrecValueType>& ls_der_samp,
                  ValueType vgs_samp,
                  ValueType weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that reduces all vector information from all processors
  /// to the root
  ///         processor
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sample_finish();

  /// Function for computing ratios of the form <f>/<g> as well as the associated
  /// variance and standard error
  void mpi_unbiased_ratio_of_means(int numSamples,
                                   std::vector<ValueType>& weights,
                                   std::vector<ValueType>& numerSamples,
                                   std::vector<ValueType>& denomSamples,
                                   ValueType& mean,
                                   ValueType& variance,
                                   ValueType& stdErr);

  /// Returns the derivatives of the cost function we are minimizing
  const std::vector<ValueType>& getAveragedDerivatives() const { return lderivs_; }

  /// helper method for updating parameter values with descent
  void updateParameters();

  /// helper method for seting step sizes for different parameter types in
  /// descent optimization
  ValueType setStepSize(int i);

  /// stores derivatives so they can be used in accelerated descent algorithm on
  /// later iterations
  void storeDerivRecord() { deriv_records_.push_back(lderivs_); }

  /// helper method for transferring information on parameter names and types to
  /// the engine
  void setupUpdate(const optimize::VariableSet& my_vars);

  /// Store a vector of parameter differences to be used by the BLM in a hybrid
  /// optimization
  void storeVectors(std::vector<ValueType>& current_params);

  /// Compute uncertainties for energy/target function and variance over a
  /// history of samples from a set of iterations
  void computeFinalizationUncertainties(std::vector<ValueType>& weights,
                                        std::vector<ValueType>& numerSamples,
                                        std::vector<ValueType>& denomSamples);

  /// Returns number of times a parameter difference vector will be stored in the
  /// optimization
  int retrieveStoreFrequency() const { return store_num_; }

  /// Returns the set of stored parameter difference vectors that will be given
  /// to the BLM
  const std::vector<std::vector<ValueType>>& retrieveHybridBLM_Input() const { return hybrid_blm_input_; }

  /// Returns the current set of parameter values
  const std::vector<ValueType>& retrieveNewParams() const { return current_params_; }

  /// Returns number of optimization steps that have been taken with descent
  int getDescentNum() const { return descent_num_; }

  /// Returns current value of omega
  ValueType getOmega() const { return omega_; }

  /// Returns value of average energy
  ValueType getEnergy() const { return e_avg_; }

  /// Returns variance of the energy
  ValueType getVariance() const { return e_var_; }

  /// Returns standard deviation of energy
  ValueType getSD() const { return e_sd_; }

  /// Returns whether an excited state is being targeted
  bool targetingExcited() const { return engine_target_excited_; }

  /// Returns the descent iteration number when on a finalizing descent section
  int getFinalDescentNum() const { return final_descent_num_; }

  /// Resets the number of vectors stored to 0 for next hybrid method
  /// macro-iteration
  void resetStorageCount() { store_count_ = 0; }

  /// Function for setting averaged derivatives, currently only used as part of a
  /// unit test of the engine's parameter update
  void setDerivs(std::vector<ValueType>& test_derivs) { lderivs_ = test_derivs; }

  /// Function for setting parameter value, used to keep descent parameter values
  /// up to date with changes that occur on BLM steps of hybrid method
  void setParamVal(int index, ValueType value) { current_params_[index] = value; }
};

} // namespace qmcplusplus
#endif
