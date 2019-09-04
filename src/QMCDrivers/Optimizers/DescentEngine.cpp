
//Prototype code for a descent engine
//Some functions modeled on those in linear method engine, but trying to avoid use of formic wrappers for things like mpi.
//Long-term(?) amibition is to eventually get rid of most of formic except for parts essential for the adaptive three shift LM and BLM.

#include<vector>
#include<string>

#include "QMCDrivers/Optimizers/DescentEngine.h"

namespace qmcplusplus
{

DescentEngine::DescentEngine(const bool targetExcited, Communicate* comm)
  : engineTargetExcited(targetExcited), myComm(comm)
{
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
 * @param q current xmlNode
 * @return true if successful
 */
/*
bool DescentEngine::parseXML(xmlNodePtr q)
{



}
*/

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
void DescentEngine::take_sample(std::vector<double> & der_rat_samp,
                                          std::vector<double> & le_der_samp,
                                          std::vector<double> & ls_der_samp,
                                          double vgs_samp,
                                          double weight_samp) {

    const size_t numOptimizables = der_rat_samp.size() - 1;
    avg_le_der_samp.resize(numOptimizables, 0.0);
    avg_der_rat_samp.resize(numOptimizables, 0.0);
    LDerivs.resize(numOptimizables, 0.0);
    
  // get the number of threads being used
 // int NumThreads = omp_get_num_threads();


  // get the thread number 
 // int myThread = omp_get_thread_num();


  e_sum += le_der_samp[0]*vgs_samp*weight_samp;

  eSquare_sum += (le_der_samp[0]*vgs_samp*weight_samp)*(le_der_samp[0]*vgs_samp*weight_samp);

  w_sum += vgs_samp*weight_samp;

  for(int i = 0; i < numOptimizables; i++)
  {
  avg_le_der_samp.at(i) += le_der_samp.at(i+1);
  avg_der_rat_samp.at(i) +=der_rat_samp.at(i+1);
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
void DescentEngine::take_sample(double local_en,
                                          double vgs_samp,
                                          double weight_samp) 
{


}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that reduces all vector information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::sample_finish() {
  
  // get rank number and number of ranks
  int my_rank = myComm->rank();
//  int num_rank = formic::mpi::size();

  // get total number of threads
 // int NumThreads = omp_get_max_threads();


    std::vector<double> etemp(3);
    
    etemp[0] = e_sum;
    etemp[1] = w_sum;
    etemp[2] = eSquare_sum;

    myComm->allreduce(etemp);
  myComm->allreduce(avg_le_der_samp);
   myComm->allreduce(avg_der_rat_samp);

   e_avg = etemp[0]/etemp[1];
   eSquare_avg = etemp[2]/etemp[1];

   w_sum = etemp[1];

   if(my_rank == 0)
   {
   std::cout << "This is e_avg: " << e_avg << std::endl;
   std::cout << "This is total weights: " << w_sum << std::endl;
   }

for (int i = 0; i < LDerivs.size();i++)
{

avg_le_der_samp.at(i) = avg_le_der_samp.at(i)/w_sum;
avg_der_rat_samp.at(i) = avg_der_rat_samp.at(i)/w_sum;

if(my_rank == 0)
{
std::cout << "Parameter # " << i << " Hamiltonian term: " << avg_le_der_samp.at(i) << std::endl;
std::cout << "Parameter # " << i <<  " Overlap term: " << avg_der_rat_samp.at(i) << std::endl;
}

if(!engineTargetExcited)
{
 LDerivs.at(i) = 2*avg_le_der_samp.at(i) - e_avg*(2*avg_der_rat_samp.at(i));
 
}



}

}

}
