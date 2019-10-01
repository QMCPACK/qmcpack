//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DESCENT_ENGINE_HEADER
#define QMCPLUSPLUS_DESCENT_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "Optimize/VariableSet.h"
#include "Configuration.h"

namespace qmcplusplus
{
class DescentEngine
{
  typedef qmcplusplus::QMCTraits::FullPrecValueType FullPrecValueType;
  typedef qmcplusplus::QMCTraits::ValueType ValueType;
  typedef qmcplusplus::QMCTraits::RealType RealType;
  typedef qmcplusplus::QMCTraits::FullPrecRealType FullPrecRealType;

private:
  //Vectors and scalars used in calculation of averaged derivatives in descent
  std::vector<FullPrecValueType> avg_le_der_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_le_der_samp_;

  std::vector<FullPrecValueType> avg_der_rat_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_der_rat_samp_;

  FullPrecValueType w_sum_;
  FullPrecValueType e_avg_;
  FullPrecValueType e_sum_;
  FullPrecValueType e_square_sum_;
  FullPrecValueType e_square_avg_;

  //Vector that stores the final averaged derivatives of the cost function
  std::vector<ValueType> lderivs_;

  //Communicator handles MPI reduction
  Communicate* my_comm_;

  //Whether to target excited state
  //Currently only ground state optimization is implemented
  bool engine_target_excited_;

  //Number of optimizable parameters
  int num_params_;


  //Vector for storing parameter values from previous optimization step
  std::vector<ValueType> params_copy_;

  //Vector for storing parameter values for current optimization step
  std::vector<ValueType> current_params_;

  //Vector for storing Lagrangian derivatives from previous optimization steps
  std::vector<std::vector<ValueType>> deriv_records_;

  //Vector for storing step size denominator values from previous optimization step
  std::vector<ValueType> denom_records_;

  //Vector for storing step size numerator values from previous optimization step
  std::vector<ValueType> numer_records_;


  //Parameter for accelerated descent recursion relation
  ValueType lambda_;
  //Vector for storing step sizes from previous optimization step.
  std::vector<ValueType> taus_;
  //Vector for storing running average of squares of the derivatives
  std::vector<ValueType> derivs_squared_;

  //Integer for keeping track of only number of descent steps taken
  int descent_num_;

  //What variety of gradient descent will be used
  std::string flavor_;

  //Step sizes for different types of parameters
  ValueType tjf_2body_eta_;
  ValueType tjf_1body_eta_;
  ValueType f_eta_;
  ValueType gauss_eta_;
  ValueType ci_eta_;
  ValueType orb_eta_;

  //Whether to gradually ramp up step sizes in descent
  bool ramp_eta_;

  //Number of steps over which to ramp up step size
  int ramp_num_;


  //Number of parameter difference vectors stored when descent is used in a hybrid optimization
  int store_num_;

  //Counter of how many vectors have been stored so far
  int store_count_;

  //Vectors of parameter names and types, used in the assignment of step sizes
  std::vector<std::string> engine_param_names_;
  std::vector<int> engine_param_types_;


  //Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<ValueType> params_for_diff_;

  //Vector for storing the input vectors to the BLM steps of hybrid method
  std::vector<std::vector<ValueType>> hybrid_blm_input_;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

public:
  //Constructor for engine
  DescentEngine(Communicate* comm, const xmlNodePtr cur);

  //Prepare for taking samples
  void prepareStorage(const int num_replicas, const int num_optimizables);

  //Sets value of average local energy
  void setEtemp(const std::vector<FullPrecRealType>& etemp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that Take Sample Data from the Host Code
  ///
  /// \param[in]  der_rat_samp   <n|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
  /// \param[in]  le_der_samp    <n|H|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
  /// \param[in]  ls_der_samp    <|S^2|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
  /// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
  /// \param[in]  weight_samp    weight for this sample
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void takeSample(const int replica_id,
                  const std::vector<FullPrecValueType>& der_rat_samp,
                  const std::vector<FullPrecValueType>& le_der_samp,
                  const std::vector<FullPrecValueType>& ls_der_samp,
                  FullPrecValueType vgs_samp,
                  FullPrecValueType weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that Take Sample Data from the Host Code
  ///
  /// \param[in]  local_en       local energy
  /// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
  /// \param[in]  weight_samp    weight for this sample
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void takeSample(FullPrecValueType local_en, FullPrecValueType vgs_samp, FullPrecValueType weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that reduces all vector information from all processors to the root
  ///         processor
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sample_finish();

  //Returns the derivatives of the cost function we are minimizing
  const std::vector<ValueType>& getAveragedDerivatives() const { return lderivs_; }

  // helper method for updating parameter values with descent
  void updateParameters();

  //helper method for seting step sizes for different parameter types in descent optimization
  ValueType setStepSize(int i);

  //stores derivatives so they can be used in accelerated descent algorithm on later iterations
  void storeDerivRecord() { deriv_records_.push_back(lderivs_); }

  //helper method for transferring information on parameter names and types to the engine
  void setupUpdate(const optimize::VariableSet& my_vars);

  //Store a vector of parameter differences to be used by the BLM in a hybrid optimization
  void storeVectors(std::vector<ValueType>& current_params);

  //Returns number of times a parameter difference vector will be stored in the optimization
  int retrieveStoreFrequency() const { return store_num_; }

  //Returns the set of stored parameter difference vectors that will be given to the BLM
  const std::vector<std::vector<ValueType>>& retrieveHybridBLM_Input() const { return hybrid_blm_input_; }

  //Returns the current set of parameter values
  const std::vector<ValueType>& retrieveNewParams() const { return current_params_; }

  //Returns number of optimization steps that have been taken with descent
  int getDescentNum() const { return descent_num_; }

  ///Function for setting averaged derivatives, currently only used as part of a unit test of the engine's parameter update
  void setDerivs(std::vector<ValueType>& test_derivs) { lderivs_ = test_derivs; }
};

} // namespace qmcplusplus
#endif
