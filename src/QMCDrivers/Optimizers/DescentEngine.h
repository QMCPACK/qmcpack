
//Prototype code for an engine to handle descent optimization


#ifndef QMCPLUSPLUS_DESCENT_ENGINE_HEADER
#define QMCPLUSPLUS_DESCENT_ENGINE_HEADER

#include <vector>
#include "Message/Communicate.h"

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

public:
  //Constructor for engine
  DescentEngine(const bool targetExcited, Communicate* comm);


  ///process xml node
  //bool processXML(xmlNodePtr cur);
  void clear_samples(const size_t numOptimizables);

  void setEtemp(std::vector<double> etemp);

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
  void take_sample(std::vector<double>& der_rat_samp,
                   std::vector<double>& le_der_samp,
                   std::vector<double>& ls_der_samp,
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
};

} // namespace qmcplusplus
#endif
