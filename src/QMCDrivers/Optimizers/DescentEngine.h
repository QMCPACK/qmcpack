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

  std::vector<FullPrecValueType> avg_numer_der_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_numer_der_samp_;

  std::vector<FullPrecValueType> avg_denom_der_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_denom_der_samp_;

  ValueType w_sum_;
  ValueType e_avg_;
  ValueType e_sum_;
  ValueType e_square_sum_;
  ValueType e_square_avg_;


  ValueType e_var_;
  ValueType e_sd_;

  //ValueType numer_sum_;
  //ValueType denom_sum_;

  ValueType numer_avg_;
  ValueType numer_var_;
  ValueType denom_avg_;
  ValueType denom_var_;
  ValueType target_avg_;
  ValueType target_var_;

  /*
  ValueType sf_; 
  ValueType sg_; 
  ValueType mf_; 
  ValueType mg_; 
  ValueType mp_; 
  ValueType ns_; 
  ValueType vf_; 
  ValueType vg_; 
  ValueType cv_; 
  ValueType other_avg_;
  ValueType other_var_;

  ValueType tsf_;
  ValueType tsg_;
  ValueType tmf_;
  ValueType tmg_;
  ValueType tmp_;
  ValueType tns_;
  ValueType tvf_;
  ValueType tvg_;
  ValueType tcv_;
  ValueType other_target_;
  ValueType other_target_var_;
*/


//Iteration to start collecting samples for final average and error blocking analysis
    int collection_step_;

      int compute_step_;

      bool collect_count_ = false;

      int final_descent_num_ = 0;

    /// \brief [in] history of sampled local energies 
    //std::vector<ValueType> _le_history;

    /// \brief [in] history of sampled |value/guiding|^2 ratios 
    std::vector<ValueType> vg_history_;
    std::vector<ValueType> final_vg_history_;

    /// \brief [in] history of sampled configuration weight 
    std::vector<ValueType> w_history_;
    std::vector<ValueType> final_w_history_;

    /// \brief a history of sampled local energies times the |value/guiding|^2 raitos
    std::vector<ValueType> lev_history_;
    
    std::vector<ValueType> final_lev_history_;

    /// \brief a history of sampled target function numerator
    //std::vector<ValueType> _tn_history;

    /// \brief a history of target function numerator times the |value/guiding|^2 ratios
    std::vector<ValueType> tnv_history_;
    std::vector<ValueType> final_tnv_history_;

    /// \brief a history of target function denominator
    //std::vector<ValueType> _td_history;

    /// \brief a history of target function denominator times the |value/guiding|^2 ratios
    std::vector<ValueType> tdv_history_;
    
    std::vector<ValueType> final_tdv_history_;

    /// \brief a history of sampled local energy square 
    //std::vector<ValueType> _les_history;

    /// \brief a history of sampled local energy square times |value/guiding|^2 ratios 
    //std::vector<ValueType> _lesv_history;



  //Vector that stores the final averaged derivatives of the cost function
  std::vector<ValueType> lderivs_;

  //Communicator handles MPI reduction
  Communicate* my_comm_;

  //Whether to target excited state
  //Currently only ground state optimization is implemented
  bool engine_target_excited_;

//Whether to use <(omega - H)^2> for targeting excited state,default is to use target above functional instead
 // bool target_excited_closest_;

 //Whether to clip samples with local energy outliers
  bool use_clipping_;

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

  //Value of omega in excited state functional
  ValueType omega_;

  //the iteration where the omega_shift parameter starts being updated
 // int update_omega_iter_;
//the number of iterations over which omega_shift is updated
 // int update_omega_steps_;

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

void setTarget(const std::vector<FullPrecRealType>& targetSums);

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
  //void takeSample(FullPrecValueType local_en, FullPrecValueType vgs_samp, FullPrecValueType weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that reduces all vector information from all processors to the root
  ///         processor
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sample_finish();


  void mpi_unbiased_ratio_of_means(int numSamples, std::vector<ValueType>& weights, std::vector<ValueType>& numerSamples,std::vector<ValueType>& denomSamples, ValueType& mean, ValueType& variance);

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

 //Adjust omega during target above functional calculation
 // void changeOmega();

  //Compute final averages for energy and variance over a history of samples from a set of iterations
  void computeFromHistory();

ValueType helperHistoryCompute(std::vector<FullPrecRealType>& weights, std::vector<FullPrecRealType>& numerator, std::vector<FullPrecRealType>& denominator,bool computing_target);

ValueType helperErrorCompute(std::vector<FullPrecRealType>& weights, std::vector<FullPrecRealType>& numerator, std::vector<FullPrecRealType>& denominator);


  //Returns number of times a parameter difference vector will be stored in the optimization
  int retrieveStoreFrequency() const { return store_num_; }

  //Returns the set of stored parameter difference vectors that will be given to the BLM
  const std::vector<std::vector<ValueType>>& retrieveHybridBLM_Input() const { return hybrid_blm_input_; }

  //Returns the current set of parameter values
  const std::vector<ValueType>& retrieveNewParams() const { return current_params_; }

  //Returns number of optimization steps that have been taken with descent
  int getDescentNum() const { return descent_num_; }

//Returns whether clipping is being used
  bool getClipping() const {return use_clipping_;}

  //Returns current value of omega
  ValueType getOmega() const {return omega_;}

  //Returns value of average energy
  ValueType getEnergy() const {return e_avg_;}

  //Returns variance of the energy
  ValueType getVariance() const {return e_var_;}

  //Returns standard deviation of energy
  ValueType getSD() const {return e_sd_;}

  //Returns whether an excited state is being targeted
  bool targetingExcited() const {return engine_target_excited_;}

  //Returns whether an adaptive omega is being used
//  bool varyingOmega() const {return update_omega_iter_ > -1;}

  int getFinalDescentNum() const {return final_descent_num_;}

  //Resets the number of vectors stored to 0 for next hybrid method macro-iteration
  void resetStorageCount() {store_count_ = 0;}

  //Function for setting averaged derivatives, currently only used as part of a unit test of the engine's parameter update
  void setDerivs(std::vector<ValueType>& test_derivs) { lderivs_ = test_derivs; }

//Function for setting parameter value, used to keep descent parameter values up to date with changes that occur on BLM steps of hybrid method
  void setParamVal(int index, ValueType value) {current_params_[index] = value;}

};

} // namespace qmcplusplus
#endif
