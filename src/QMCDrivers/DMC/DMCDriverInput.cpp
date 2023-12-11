//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "DMCDriverInput.h"

namespace qmcplusplus
{
DMCDriverInput::DMCDriverInput(int walkers_per_rank) {}

void DMCDriverInput::readXML(xmlNodePtr node)
{
  ParameterSet parameter_set_;
  std::string reconfig_str;
  std::string refE_update_scheme_str;
  parameter_set_.add(reconfig_str, "reconfiguration", {"no", "yes", "runwhileincorrect"});
  parameter_set_.add(NonLocalMove, "nonlocalmove", {"no", "yes", "v0", "v1", "v3"});
  parameter_set_.add(NonLocalMove, "nonlocalmoves", {"no", "yes", "v0", "v1", "v3"});
  parameter_set_.add(max_age_, "MaxAge");
  parameter_set_.add(feedback_, "feedback");
  parameter_set_.add(refE_update_scheme_str, "refenergy_update_scheme", {"unlimited_history", "limited_history"});

  // from DMC.cpp put(xmlNodePtr)
  parameter_set_.add(branch_interval_, "branchInterval");
  parameter_set_.add(branch_interval_, "branchinterval");
  parameter_set_.add(branch_interval_, "substeps");
  parameter_set_.add(branch_interval_, "subStep");
  parameter_set_.add(branch_interval_, "sub_stepd");

  //from NonLocalTOperator.cpp
  parameter_set_.add(alpha_, "alpha");
  parameter_set_.add(gamma_, "gamma");

  parameter_set_.add(reserve_, "reserve");

  parameter_set_.put(node);

  if (reconfig_str == "yes")
    throw std::runtime_error("Reconfiguration is currently broken and gives incorrect results. Use dynamic "
                             "population control by setting reconfiguration=\"no\" or removing the reconfiguration "
                             "option from the DMC input section. If accessing the broken reconfiguration code path "
                             "is still desired, set reconfiguration to \"runwhileincorrect\" instead of \"yes\".");
  reconfiguration_ = (reconfig_str == "yes");

  if (NonLocalMove == "yes" || NonLocalMove == "v0")
    app_summary() << "  Using Non-local T-moves v0, M. Casula, PRB 74, 161102(R) (2006)";
  else if (NonLocalMove == "v1")
    app_summary() << "  Using Non-local T-moves v1, M. Casula et al., JCP 132, 154113 (2010)";
  else if (NonLocalMove == "v3")
    app_summary() << "  Using Non-local T-moves v3, an approximation to v1";
  else
    app_summary() << "  Using Locality Approximation";
  app_summary() << std::endl;

  // TODO: similar check for alpha and gamma
  if (max_age_ < 0)
    throw std::runtime_error("Illegal input for MaxAge in DMC input section");
  if (branch_interval_ < 0)
    throw std::runtime_error("Illegal input for branchInterval or substeps in DMC input section");

  if (reserve_ < 1.0)
    throw std::runtime_error("You can only reserve walkers above the target walker count");

  if (refE_update_scheme_str == "unlimited_history")
    refenergy_update_scheme_ = DMCRefEnergyScheme::UNLIMITED_HISTORY;
  else
    refenergy_update_scheme_ = DMCRefEnergyScheme::LIMITED_HISTORY;
}

std::ostream& operator<<(std::ostream& o_stream, const DMCDriverInput& dmci) { return o_stream; }

} // namespace qmcplusplus
