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


#include "WFOptDriverInput.h"
#include "io/OhmmsData/libxmldefs.h"


namespace qmcplusplus
{
void WFOptDriverInput::readXML(xmlNodePtr xml_input)
{
  ParameterSet parameter_set;
  parameter_set.add(opt_crowd_size_, "opt_crowd_size");
  parameter_set.add(opt_num_crowds_, "opt_num_crowds");
  parameter_set.put(xml_input);

  opt_xml_node_ = xml_input->children;

  xmlNodePtr cur = xml_input->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname.find("optimize") < cname.size())
    {
      const std::string att(getXMLAttributeValue(cur, "method"));
      if (!att.empty())
        opt_method_ = att;
      opt_xml_node_ = cur;
    }
    cur = cur->next;
  }
}

} // namespace qmcplusplus
