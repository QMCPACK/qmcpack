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

#ifndef QMCPLUSPLUS_DMCDRIVERINPUT_H
#define QMCPLUSPLUS_DMCDRIVERINPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
/** Input representation for DMC driver class runtime parameters
 */
class DMCDriverInput
{
public:
  using IndexType             = QMCTraits::IndexType;
  using RealType              = QMCTraits::RealType;
  using FullPrecisionRealType = QMCTraits::FullPrecRealType;
  DMCDriverInput(){};
  DMCDriverInput(int walkers_per_rank);
  void readXML(xmlNodePtr xml_input);

  bool get_reconfiguration() const { return reconfiguration_; }
  IndexType get_max_age() const { return max_age_; }
  IndexType get_branch_interval() const { return branch_interval_; }

private:
  /** @ingroup Parameters for DMC Driver
   *  @{
   *  
   *  Do not write out blocks of gets for variables like this
   *  there is are code_generation tools in QMCPACK_ROOT/utils/code_tools
   */
  ///Interval between branching
  IndexType branch_interval_ = 1;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  /// reconfiguration flag
  bool reconfiguration_ = true;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  ///input std::string to use fast gradient
  std::string UseFastGrad;
  ///input to control maximum age allowed for walkers.
  IndexType max_age_ = 0;
public:

  friend std::ostream& operator<<(std::ostream& o_stream, const DMCDriverInput& vmci);
};

extern std::ostream& operator<<(std::ostream& o_stream, const DMCDriverInput& vmci);

} // namespace qmcplusplus
#endif
