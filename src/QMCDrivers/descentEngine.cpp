
//Prototype code for a descent engine
//Some functions modeled on those in linear method engine, but trying to avoid use of formic wrappers for things like mpi.
//Long-term(?) amibition is to eventually get rid of most of formic except for parts essential for the adaptive three shift LM and BLM.

#include<vector>
#include<string>

#include "descentEngine.h"

//#include<formic/utils/openmp.h>


//#include "Message/MPIObjectBase.h"
//#include "Message/CommOperators.h"






cqmc::engine::descentEngine::descentEngine(const int numParams, const bool targetExcited)
{

    numOptimizables = numParams;
    avg_le_der_samp.resize(numOptimizables, 0.0);
    avg_der_rat_samp.resize(numOptimizables, 0.0);

    LDerivs.resize(numOptimizables, 0.0);

    engineTargetExcited = targetExcited;


}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
 * @param q current xmlNode
 * @return true if successful
 */
/*
bool cqmc::engine::descentEngine::parseXML(xmlNodePtr q)
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
void cqmc::engine::descentEngine::take_sample(std::vector<double> & der_rat_samp,
                                          std::vector<double> & le_der_samp,
                                          std::vector<double> & ls_der_samp,
                                          double vgs_samp,
                                          double weight_samp) {

    
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
void cqmc::engine::descentEngine::take_sample(double local_en,
                                          double vgs_samp,
                                          double weight_samp) 
{


}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that reduces all vector information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::descentEngine::sample_finish() {
  
  // get rank number and number of ranks
//  int my_rank = formic::mpi::rank();
//  int num_rank = formic::mpi::size();

  // get total number of threads
 // int NumThreads = omp_get_max_threads();


  myComm->allreduce(avg_le_der_samp);
   myComm->allreduce(avg_der_rat_samp);

   e_avg = e_sum/w_sum;

for (int i = 0; i < LDerivs.size();i++)
{

avg_le_der_samp.at(i) = avg_le_der_samp.at(i)/w_sum;
avg_der_rat_samp.at(i) = avg_der_rat_samp.at(i)/w_sum;

if(!engineTargetExcited)
{
 LDerivs.at(i) = 2*avg_le_der_samp.at(i) - e_avg*(2*avg_der_rat_samp.at(i));
}



}

}
