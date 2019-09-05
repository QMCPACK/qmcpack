
//Prototype code for a descent engine
//Some functions modeled on those in linear method engine, but trying to avoid use of formic wrappers for things like mpi.
//Long-term(?) amibition is to eventually get rid of most of formic except for parts essential for the adaptive three shift LM and BLM.

#include <vector>
#include <string>

#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{
DescentEngine::DescentEngine(const bool targetExcited, Communicate* comm)
    : engineTargetExcited(targetExcited), myComm(comm)
{}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
 * @param q current xmlNode
 * @return true if successful
 */
/*
bool DescentEngine::parseXML(xmlNodePtr q)
{



}
*/

void DescentEngine::clear_samples(const size_t numOptimizables)
{
  avg_le_der_samp.resize(numOptimizables);
  avg_der_rat_samp.resize(numOptimizables);
  LDerivs.resize(numOptimizables);

  std::fill(avg_le_der_samp.begin(), avg_le_der_samp.end(), 0.0);
  std::fill(avg_der_rat_samp.begin(), avg_der_rat_samp.end(), 0.0);

  w_sum       = 0;
  e_avg       = 0;
  e_sum       = 0;
  eSquare_sum = 0;
  eSquare_avg = 0;
}

void DescentEngine::setEtemp(std::vector<double> etemp)
{
    e_sum = etemp[0];
    w_sum = etemp[1];
    eSquare_sum = etemp[2];
    e_avg = e_sum/w_sum;
    eSquare_avg = eSquare_sum/w_sum;

    int my_rank = myComm->rank();
    if(my_rank == 0)
    {
    std::cout << "e_sum: " << e_sum << std::endl;
    std::cout << "w_sum: " << w_sum << std::endl;
    std::cout << "e_avg: " << e_avg << std::endl;
    std::cout << "eSquare_sum: " << eSquare_sum << std::endl;
    std::cout << "eSquare_avg: " << eSquare_avg << std::endl;
    }
}

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
void DescentEngine::take_sample(std::vector<double>& der_rat_samp,
                                std::vector<double>& le_der_samp,
                                std::vector<double>& ls_der_samp,
                                double vgs_samp,
                                double weight_samp)
{
  const size_t numOptimizables = der_rat_samp.size() - 1;

  


  for (int i = 0; i < numOptimizables; i++)
  {
    avg_le_der_samp.at(i) += le_der_samp.at(i + 1);
    avg_der_rat_samp.at(i) += der_rat_samp.at(i + 1);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
///
/// \param[in]  local_en       local energy
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::take_sample(double local_en, double vgs_samp, double weight_samp) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that reduces all vector information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::sample_finish()
{
  // get rank number and number of ranks
  int my_rank = myComm->rank();

  myComm->allreduce(avg_le_der_samp);
  myComm->allreduce(avg_der_rat_samp);

  for (int i = 0; i < LDerivs.size(); i++)
  {
    avg_le_der_samp.at(i)  = avg_le_der_samp.at(i) / w_sum;
    avg_der_rat_samp.at(i) = avg_der_rat_samp.at(i) / w_sum;

    if (my_rank == 0)
    {
      std::cout << "Parameter # " << i << " Hamiltonian term: " << avg_le_der_samp.at(i) << std::endl;
      std::cout << "Parameter # " << i << " Overlap term: " << avg_der_rat_samp.at(i) << std::endl;
    }

    if (!engineTargetExcited)
    {
      LDerivs.at(i) = 2 * avg_le_der_samp.at(i) - e_avg * (2 * avg_der_rat_samp.at(i));
    }
  }
}

} // namespace qmcplusplus
