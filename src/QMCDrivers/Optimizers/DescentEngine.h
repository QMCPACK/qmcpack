
//Prototype code for an engine to handle descent optimization


#ifndef QMCPLUSPLUS_DESCENT_ENGINE_HEADER
#define QMCPLUSPLUS_DESCENT_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "Optimize/VariableSet.h"

namespace qmcplusplus
{
class DescentEngine
{
private:
  std::vector<double> avg_le_der_samp;
  std::vector<double> avg_der_rat_samp;

  double w_sum;
  double e_avg;
  double e_sum;
  double eSquare_sum;
  double eSquare_avg;

  Communicate* myComm;

  std::vector<double> LDerivs;

  bool engineTargetExcited;

  int numParams;


  //Vector for storing parameter values from previous optimization step
  std::vector<double> paramsCopy;

  //Vector for storing parameter values for current optimization step
  std::vector<double> currentParams;

  //Vector for storing Lagrangian derivatives from previous optimization steps
  std::vector<std::vector<double>> derivRecords;

  //Vector for storing step size denominator values from previous optimization step
  std::vector<double> denomRecords;

  //Vector for storing step size numerator values from previous optimization step
  std::vector<double> numerRecords;


  //Parameter for accelerated descent recursion relation
  double lambda;
  //Vector for storing step sizes from previous optimization step.
  std::vector<double> taus;
  //Vector for storing running average of squares of the derivatives
  std::vector<double> derivsSquared;
  //Integer for keeping track of the iteration number
  int stepNum;
  //Integer for keeping track of only number of descent steps taken
  int descentNum;

  //What variety of gradient descent will be used
  std::string flavor;

  //Step sizes for different types of parameters
  double TJF_2Body_eta;
  double TJF_1Body_eta;
  double F_eta;
  double Gauss_eta;
  double CI_eta;
  double Orb_eta;

  //Whether to gradually ramp up step sizes in descent
  bool ramp_eta;

  //Number of steps over which to ramp up step size
  int ramp_num;


  std::vector<std::string> engineParamNames;

  std::vector<int> engineParamTypes;


  //Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<double> paramsForDiff;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

public:
  //Constructor for engine
  DescentEngine(Communicate* comm, const xmlNodePtr cur);

  void clear_samples(const int numOptimizables);

  void setEtemp(const std::vector<double>& etemp);

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
  void take_sample(const std::vector<double>& der_rat_samp,
                   const std::vector<double>& le_der_samp,
                   const std::vector<double>& ls_der_samp,
                   double vgs_samp,
                   double weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that Take Sample Data from the Host Code
  ///
  /// \param[in]  local_en       local energy
  /// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
  /// \param[in]  weight_samp    weight for this sample
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void take_sample(double local_en, double vgs_samp, double weight_samp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Function that reduces all vector information from all processors to the root
  ///         processor
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  void sample_finish();

  const std::vector<double>& getAveragedDerivatives() { return LDerivs; }

  // helper method for updating parameter values with descent
  void updateParameters(int stepNum, int descentNum);

  //helper method for seting step sizes for different parameter types in descent optimization
  double setStepSize(int i);

  void storeDerivRecord() { derivRecords.push_back(LDerivs); }

  void setupUpdate(const optimize::VariableSet& myVars);

  std::vector<double> retrieveNewParams() { return currentParams; }
};

} // namespace qmcplusplus
#endif
