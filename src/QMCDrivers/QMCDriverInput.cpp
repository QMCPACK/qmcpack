//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  
 */

#include "OhmmsData/AttributeSet.h"
#include "QMCDriverInput.h"
#include "Concurrency/Info.hpp"

namespace qmcplusplus
{
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
  // ParameterSet has an dependency on the lifetime of the backing xmlNodePtr
  // so its better it not live long
  ParameterSet parameter_set;
  parameter_set.add(RollBackBlocks_, "rewind", "int");
  parameter_set.add(store_config_period_, "storeconfigs", "int");
  parameter_set.add(store_config_period_, "store_configs", "int");
  parameter_set.add(recalculate_properties_period_, "checkProperties", "int");
  parameter_set.add(recalculate_properties_period_, "checkproperties", "int");
  parameter_set.add(recalculate_properties_period_, "check_properties", "int");
  parameter_set.add(config_dump_period_.period, "recordconfigs", "int");
  parameter_set.add(config_dump_period_.period, "record_configs", "int");
  parameter_set.add(starting_step_, "current", "int");
  parameter_set.add(max_blocks_, "blocks", "int");
  parameter_set.add(max_steps_, "steps", "int");
  parameter_set.add(sub_steps_, "substeps", "int");
  parameter_set.add(sub_steps_, "sub_steps", "int");
  parameter_set.add(warmup_steps_, "warmupsteps", "int");
  parameter_set.add(warmup_steps_, "warmup_steps", "int");
  parameter_set.add(num_crowds_, "crowds", "int");
  parameter_set.add(walkers_per_rank_, "walkers_per_rank", "int");
  parameter_set.add(total_walkers_, "total_walkers", "int");
  parameter_set.add(steps_between_samples_, "stepsbetweensamples", "int");
  parameter_set.add(samples_per_thread_, "samplesperthread", "real");
  parameter_set.add(requested_samples_, "samples", "real");
  parameter_set.add(tau_, "timestep", "AU");
  parameter_set.add(tau_, "time_step", "AU");
  parameter_set.add(tau_, "tau", "AU");
  parameter_set.add(max_cpu_secs_, "maxcpusecs", "real");
  parameter_set.add(blocks_between_recompute_, "blocks_between_recompute", "int");
  parameter_set.add(drift_modifier_, "drift_modifier", "string");
  parameter_set.add(drift_modifier_unr_a_, "drift_UNR_a", "double");
  parameter_set.add(max_disp_sq_, "maxDisplSq", "double");

  OhmmsAttributeSet aAttrib;
  // first stage in from QMCDriverFactory
  aAttrib.add(qmc_method_, "method");
  aAttrib.add(update_mode_, "move");
  aAttrib.add(scoped_profiling_, "profiling");

  aAttrib.add(k_delay_, "kdelay");
  // This does all the parameter parsing setup in the constructor
  aAttrib.put(cur);
  if (cur != NULL)
  {
    //initialize the parameter set
    parameter_set.put(cur);

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
        dump_config_ = (check_point_period_.period > 0);
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

} // namespace qmcplusplus
