//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: XMLParsingString.hpp
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 */
#include "XMLParsingString.h"
#include "libxmldefs.h"

std::string getXMLAttributeValue(const xmlNodePtr cur, const std::string_view name)
{
  std::string attr_value;
  xmlChar* attr_char = xmlGetProp(cur, castCharToXMLChar(name.data()));
  if (attr_char)
  {
    attr_value.assign(castXMLCharToChar(attr_char));
    xmlFree(attr_char);
  }
  return attr_value;
}
