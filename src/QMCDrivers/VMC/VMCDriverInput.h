//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VMCDRIVERINPUT_H
#define QMCPLUSPLUS_VMCDRIVERINPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
class VMCDriverInput
{
public:
  using IndexType             = QMCTraits::IndexType;
  using RealType              = QMCTraits::RealType;
  using FullPrecisionRealType = QMCTraits::FullPrecRealType;

  VMCDriverInput(){};
  VMCDriverInput(int walkers_per_rank, const std::string& use_drift);
  void readXML(xmlNodePtr& xml_input);

protected:
  ///store any parameter that has to be read from a file
  //ParameterSet parameter_set_;

  /** @ingroup Parameters for VMC Driver
   *  @{
   *  All unshared input should be here
   *  
   *  Do not write out blocks of gets for variables like this
   *  there is are code_generation tools in QMCPACK_ROOT/utils/code_generation
   */

  IndexType requested_walkers_per_rank_ = 0;
  std::string use_drift_{"yes"};

  /** @} */

public:
  IndexType get_requested_walkers_per_rank() const { return requested_walkers_per_rank_; }
  const std::string get_use_drift() const { return use_drift_; }
};

} // namespace qmcplusplus
#endif
