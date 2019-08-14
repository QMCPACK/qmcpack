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

#include "QMCDrivers/QMCDriverInput.h"

namespace qmcplusplus
{
class VMCDriverInput
{
public:
  VMCDriverInput(int qmc_section_count) : qmcdriver_input(qmc_section_count) {}
  inline QMCDriverInput& get_qmcdriver_input() { return qmcdriver_input; }
  void readXML(xmlNodePtr& xml_input);
protected:
  QMCDriverInput qmcdriver_input;
};

} // namespace qmcplusplus
#endif
