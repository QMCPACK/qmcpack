//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WFOPTDRIVERINPUT_H
#define QMCPLUSPLUS_WFOPTDRIVERINPUT_H

#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
class WFOptDriverInput
{
public:
  WFOptDriverInput() {}
  void readXML(xmlNodePtr xml_input);

protected:
  int opt_crowd_size_ = 1;
  int opt_num_crowds_ = 0;

  std::string opt_method_;

  xmlNodePtr opt_xml_node_;

public:
  int get_opt_crowd_size() const { return opt_crowd_size_; }
  int get_opt_num_crowds() const { return opt_num_crowds_; }

  const std::string& get_opt_method() const { return opt_method_; }

  xmlNodePtr get_opt_xml_node() const { return opt_xml_node_; }
};


} // namespace qmcplusplus

#endif
