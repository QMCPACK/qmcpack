//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "libxmldefs.h"
#include "ModernStringUtils.hpp"

std::string getNodeName(xmlNodePtr cur)
{
  return qmcplusplus::lowerCase(castXMLCharToChar(cur->name));
}
