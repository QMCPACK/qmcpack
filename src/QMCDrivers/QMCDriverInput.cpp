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

  std::string serialize_walkers;
  std::string debug_checks_str;

  ParameterSet parameter_set;
  parameter_set.add(store_config_period_, "storeconfigs");
  parameter_set.add(store_config_period_, "store_configs");
  parameter_set.add(recalculate_properties_period_, "checkProperties");
  parameter_set.add(recalculate_properties_period_, "checkproperties");
  parameter_set.add(recalculate_properties_period_, "check_properties");
  parameter_set.add(config_dump_period_.period, "recordconfigs");
  parameter_set.add(config_dump_period_.period, "record_configs");
  parameter_set.add(starting_step_, "current");
  parameter_set.add(max_blocks_, "blocks");
  parameter_set.add(max_steps_, "steps");
  parameter_set.add(sub_steps_, "substeps");
  parameter_set.add(sub_steps_, "sub_steps");
  parameter_set.add(warmup_steps_, "warmupsteps");
  parameter_set.add(warmup_steps_, "warmup_steps");
  parameter_set.add(num_crowds_, "crowds");
  parameter_set.add(serialize_walkers, "crowd_serialize_walkers", {"no", "yes"});
  parameter_set.add(walkers_per_rank_, "walkers_per_rank");
  parameter_set.add(walkers_per_rank_, "walkers", {}, TagStatus::UNSUPPORTED);
  parameter_set.add(total_walkers_, "total_walkers");
  parameter_set.add(steps_between_samples_, "stepsbetweensamples", {}, TagStatus::UNSUPPORTED);
  parameter_set.add(samples_per_thread_, "samplesperthread", {}, TagStatus::UNSUPPORTED);
  parameter_set.add(requested_samples_, "samples");
  parameter_set.add(tau_, "timestep");
  parameter_set.add(tau_, "time_step");
  parameter_set.add(tau_, "tau");
  parameter_set.add(spin_mass_, "spin_mass");
  parameter_set.add(blocks_between_recompute_, "blocks_between_recompute");
  parameter_set.add(drift_modifier_, "drift_modifier");
  parameter_set.add(drift_modifier_unr_a_, "drift_UNR_a");
  parameter_set.add(max_disp_sq_, "maxDisplSq");
  parameter_set.add(debug_checks_str, "debug_checks",
                    {"no", "all", "checkGL_after_load", "checkGL_after_moves", "checkGL_after_tmove"});

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

  crowd_serialize_walkers_ = serialize_walkers == "yes";
  if (crowd_serialize_walkers_)
    app_summary() << "  Batched operations are serialized over walkers." << std::endl;
  if (scoped_profiling_)
    app_summary() << "  Profiler data collection is enabled in this driver scope." << std::endl;

  if (debug_checks_str == "no")
    debug_checks_ = DriverDebugChecks::ALL_OFF;
  else
  {
    if (debug_checks_str == "all" || debug_checks_str == "checkGL_after_load")
      debug_checks_ |= DriverDebugChecks::CHECKGL_AFTER_LOAD;
    if (debug_checks_str == "all" || debug_checks_str == "checkGL_after_moves")
      debug_checks_ |= DriverDebugChecks::CHECKGL_AFTER_MOVES;
    if (debug_checks_str == "all" || debug_checks_str == "checkGL_after_tmove")
      debug_checks_ |= DriverDebugChecks::CHECKGL_AFTER_TMOVE;
  }

  if (check_point_period_.period < 1)
    check_point_period_.period = max_blocks_;
}

} // namespace qmcplusplus
