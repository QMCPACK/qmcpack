//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  
 */

#include "OhmmsData/AttributeSet.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "Concurrency/Info.hpp"

namespace qmcplusplus
{
QMCDriverInput::QMCDriverInput(int qmc_section_count) : qmc_section_count_(qmc_section_count)
{
  parameter_set_.add(RollBackBlocks_, "rewind", "int");
  parameter_set_.add(store_config_period_, "storeconfigs", "int");
  parameter_set_.add(store_config_period_, "store_configs", "int");
  parameter_set_.add(recalculate_properties_period_, "checkProperties", "int");
  parameter_set_.add(recalculate_properties_period_, "checkproperties", "int");
  parameter_set_.add(recalculate_properties_period_, "check_properties", "int");
  parameter_set_.add(config_dump_period_.period, "recordconfigs", "int");
  parameter_set_.add(config_dump_period_.period, "record_configs", "int");
  parameter_set_.add(starting_step_, "current", "int");
  parameter_set_.add(max_blocks_, "blocks", "int");
  parameter_set_.add(max_steps_, "steps", "int");
  parameter_set_.add(sub_steps_, "substeps", "int");
  parameter_set_.add(sub_steps_, "sub_steps", "int");
  parameter_set_.add(warmup_steps_, "warmupsteps", "int");
  parameter_set_.add(warmup_steps_, "warmup_steps", "int");
  parameter_set_.add(requested_walkers_per_rank_, "walkers", "int");
  parameter_set_.add(num_crowds_, "crowds", "int");
  parameter_set_.add(steps_between_samples_, "stepsbetweensamples", "int");
  parameter_set_.add(samples_per_thread_, "samplesperthread", "real");
  parameter_set_.add(samples_per_thread_, "dmcwalkersperthread", "real");
  parameter_set_.add(requested_samples_, "samples", "real");
  parameter_set_.add(tau_, "timestep", "AU");
  parameter_set_.add(tau_, "time_step", "AU");
  parameter_set_.add(tau_, "tau", "AU");
  parameter_set_.add(max_cpu_secs_, "maxcpusecs", "real");
  parameter_set_.add(blocks_between_recompute_, "blocks_between_recompute", "int");
}

/** Reads qmc section xml node parameters
 *
 * All shared parameters are read here
 * attribute list
 * - checkpoint="-1|0|n" default=-1
 *   -- 1 = do not write anything
 *   -- 0 = dump after the completion of a qmc section
 *   -- n = dump after n blocks
 * - kdelay = "0|1|n" default=0
 */
void QMCDriverInput::readXML(xmlNodePtr cur)
{
  OhmmsAttributeSet aAttrib;

  // first stage in from QMCDriverFactory
  aAttrib.add(qmc_method_, "method");
  aAttrib.add(update_mode_, "move");


  aAttrib.add(k_delay_, "kdelay");
  // This does all the parameter parsing setup in the constructor
  aAttrib.put(cur);
  if (cur != NULL)
  {
    //initialize the parameter set
    parameter_set_.put(cur);

    xmlNodePtr tcur = cur->children;
    //determine how often to print walkers to hdf5 file
    while (tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if (cname == "record")
      {
        //dump walkers for optimization
        OhmmsAttributeSet rAttrib;
        rAttrib.add(walker_dump_period_.stride, "stride");
        rAttrib.add(walker_dump_period_.period, "period");
        rAttrib.put(tcur);
      }
      else if (cname == "checkpoint")
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(check_point_period_.stride, "stride");
        rAttrib.add(check_point_period_.period, "period");
        rAttrib.put(tcur);
        //DumpConfig=(Period4CheckPoint>0);
      }
      else if (cname == "dumpconfig")
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(config_dump_period_.stride, "stride");
        rAttrib.add(config_dump_period_.period, "period");
        rAttrib.put(tcur);
      }
      else if (cname == "random")
      {
        reset_random_ = true;
      }
      tcur = tcur->next;
    }
  }

  if (check_point_period_.period < 1)
    check_point_period_.period = max_blocks_;
}


/** Input representation for Driver base class runtime parameters
 */
// void QMCDriverInput::readXML(xmlNodePtr cur) {
//       //grep minimumTargetWalker
//   int target_min = -1;
//   ParameterSet p;
//   p.add(target_min, "minimumtargetwalkers", "int"); //p.add(target_min,"minimumTargetWalkers","int");
//   p.add(target_min, "minimumsamples", "int");       //p.add(target_min,"minimumSamples","int");
//   p.put(q);

//   app_log() << "\n<vmc function=\"put\">"
//             << "\n  qmc_counter=" << qmc_common.qmc_counter << "  my_counter=" << MyCounter << std::endl;
//   if (qmc_common.qmc_counter && MyCounter)
//   {
//     nSteps               = prevSteps;
//     nStepsBetweenSamples = prevStepsBetweenSamples;
//   }
//   else
//   {
//     int nw = W.getActiveWalkers();
//     //compute samples and overwrite steps for the given samples
//     int Nthreads = omp_get_max_threads();
//     int Nprocs   = myComm->size();
//     //target samples set by samples or samplesperthread/dmcwalkersperthread
//     nTargetPopulation = std::max(nTargetPopulation, nSamplesPerThread * Nprocs * Nthreads);
//     nTargetSamples    = static_cast<int>(std::ceil(nTargetPopulation));

//     if (nTargetSamples)
//     {
//       int nwtot      = nw * Nprocs; //total number of walkers used by this qmcsection
//       nTargetSamples = std::max(nwtot, nTargetSamples);
//       if (target_min > 0)
//       {
//         nTargetSamples    = std::max(nTargetSamples, target_min);
//         nTargetPopulation = std::max(nTargetPopulation, static_cast<RealType>(target_min));
//       }
//       nTargetSamples = ((nTargetSamples + nwtot - 1) / nwtot) *
//           nwtot; // nTargetSamples are always multiples of total number of walkers
//       nSamplesPerThread = nTargetSamples / Nprocs / Nthreads;
//       int ns_target     = nTargetSamples * nStepsBetweenSamples; //total samples to generate
//       int ns_per_step   = Nprocs * nw;                           //total samples per step
//       nSteps            = std::max(nSteps, (ns_target / ns_per_step + nBlocks - 1) / nBlocks);
//       Period4WalkerDump = nStepsBetweenSamples = (ns_per_step * nSteps * nBlocks) / nTargetSamples;
//     }
//     else
//     {
//       Period4WalkerDump = nStepsBetweenSamples = (nBlocks + 1) * nSteps; //some positive number, not used
//       nSamplesPerThread                        = 0;
//     }
//   }
//   prevSteps               = nSteps;
//   prevStepsBetweenSamples = nStepsBetweenSamples;

//   app_log() << "  time step      = " << Tau << std::endl;
//   app_log() << "  blocks         = " << nBlocks << std::endl;
//   app_log() << "  steps          = " << nSteps << std::endl;
//   app_log() << "  substeps       = " << nSubSteps << std::endl;
//   app_log() << "  current        = " << CurrentStep << std::endl;
//   app_log() << "  target samples = " << nTargetPopulation << std::endl;
//   app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << std::endl << std::endl;
//   app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << std::endl;

//   m_param.get(app_log());

//   if (DumpConfig)
//   {
//     app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint
//               << " blocks" << std::endl;
//   }
//   else
//   {
//     app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << std::endl;
//   }

//   if (Period4WalkerDump > 0)
//     app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << std::endl;

//   app_log() << "</vmc>" << std::endl;
//   app_log().flush();

//   return true;

//}

} // namespace qmcplusplus
