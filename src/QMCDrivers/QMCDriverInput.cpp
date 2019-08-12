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

/** Parses the xml input file for parameter definitions for a single qmc simulation.
 *
 * Basic parameters are handled here and each driver will perform its own initialization with the input
 * attribute list
 * - checkpoint="-1|0|n" default=-1
 *   -- 1 = do not write anything
 *   -- 0 = dump after the completion of a qmc section
 *   -- n = dump after n blocks
 * - kdelay = "0|1|n" default=0
 */
void QMCDriverInput::putQMCInfo(xmlNodePtr cur)
{
  OhmmsAttributeSet aAttrib;
  aAttrib.add(k_delay_, "kdelay");
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
void QMCDriverInput::put(xmlNodePtr cur) {}

} // namespace qmcplusplus
