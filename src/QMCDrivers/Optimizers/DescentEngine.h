
//Prototype code for an engine to handle descent optimization


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
private:
    //Vectors and scalars used in calculation of averaged derivatives in descent
  std::vector<FullPrecValueType> avg_le_der_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_le_der_samp_;

  std::vector<FullPrecValueType> avg_der_rat_samp_;
  std::vector<std::vector<FullPrecValueType>> replica_der_rat_samp_;

  FullPrecValueType w_sum;
  FullPrecValueType e_avg;
  FullPrecValueType e_sum;
  FullPrecValueType eSquare_sum;
  FullPrecValueType eSquare_avg;

  std::vector<FullPrecValueType> LDerivs;

  //Communicator handles MPI reduction
  Communicate* myComm;

  //Whether to target excited state
  //Currently only ground state optimization is implemented
  bool engineTargetExcited;

  //Number of optimizable parameters
  int numParams;


  //Vector for storing parameter values from previous optimization step
  std::vector<ValueType> paramsCopy;

  //Vector for storing parameter values for current optimization step
  std::vector<ValueType> currentParams;

  //Vector for storing Lagrangian derivatives from previous optimization steps
  std::vector<std::vector<ValueType>> derivRecords;

  //Vector for storing step size denominator values from previous optimization step
  std::vector<ValueType> denomRecords;

  //Vector for storing step size numerator values from previous optimization step
  std::vector<ValueType> numerRecords;


  //Parameter for accelerated descent recursion relation
  ValueType lambda;
  //Vector for storing step sizes from previous optimization step.
  std::vector<ValueType> taus;
  //Vector for storing running average of squares of the derivatives
  std::vector<ValueType> derivsSquared;

  //Integer for keeping track of only number of descent steps taken
  int descent_num_;

  //What variety of gradient descent will be used
  std::string flavor;

  //Step sizes for different types of parameters
  ValueType TJF_2Body_eta;
  ValueType TJF_1Body_eta;
  ValueType F_eta;
  ValueType Gauss_eta;
  ValueType CI_eta;
  ValueType Orb_eta;

  //Whether to gradually ramp up step sizes in descent
  bool ramp_eta;

  //Number of steps over which to ramp up step size
  int ramp_num;


  //Number of parameter difference vectors stored when descent is used in a hybrid optimization
  int store_num;

  //Counter of how many vectors have been stored so far
  int store_count;

  //Vectors of parameter names and types, used in the assignment of step sizes
  std::vector<std::string> engineParamNames;
  std::vector<int> engineParamTypes;


  //Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<ValueType> paramsForDiff;

  //Vector for storing the input vectors to the BLM steps of hybrid method
  std::vector<std::vector<ValueType>> hybridBLM_Input;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

public:
  //Constructor for engine
  DescentEngine(Communicate* comm, const xmlNodePtr cur);

  void prepareStorage(const int num_replicas, const int num_optimizables);

  void setEtemp(const std::vector<FullPrecValueType>& etemp);

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

  const std::vector<ValueType>& getAveragedDerivatives() const { return LDerivs; }

  // helper method for updating parameter values with descent
  void updateParameters();

  //helper method for seting step sizes for different parameter types in descent optimization
  ValueType setStepSize(int i);

  //stores derivatives so they can be used in accelerated descent algorithm on later iterations
  void storeDerivRecord() { derivRecords.push_back(LDerivs); }

  //helper method for transferring information on parameter names and types to the engine
  void setupUpdate(const optimize::VariableSet& myVars);

  void storeVectors(std::vector<ValueType>& currentParams);

  int retrieveStoreFrequency() const { return store_num; }

  const std::vector<std::vector<ValueType>>& retrieveHybridBLM_Input() const { return hybridBLM_Input; }

  const std::vector<ValueType>& retrieveNewParams() const { return currentParams; }

  int getDescentNum() const { return descent_num_; }
};

} // namespace qmcplusplus
#endif
