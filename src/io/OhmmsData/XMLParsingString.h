//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file XMLParsingString.h
 *
 * XMLNodeString convert xmlNode contents into a std::string
 * XMLAttrString convert one xmlNode attribute into a std::string
 */
#ifndef QMCPLUSPLUS_XMLSTRING_H
#define QMCPLUSPLUS_XMLSTRING_H

#include <string>
#include <string_view>
#include <libxml/xmlmemory.h>
#include <libxml/tree.h>

/** convert xmlNode contents into a std::string
 *  \todo this should just be a function that takes a cur
 *  and returns a std::string.
 *  setXMLNodeContent is OOP without reason.
 */
class XMLNodeString : public std::string
{
public:
  /// construct a string from an xmlNode
  XMLNodeString(const xmlNodePtr cur)
  {
    xmlChar* node_char = xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1);
    if (node_char)
    {
      assign((const char*)node_char);
      xmlFree(node_char);
    }
  }

  /// expose base class constructors
  XMLNodeString(const std::string& in) : std::string(in) {}

  XMLNodeString(const char* in) : std::string(in) {}

  /// write a string to an xmlNode
  void setXMLNodeContent(xmlNodePtr cur) const { xmlNodeSetContent(cur, (const xmlChar*)(this->c_str())); }
};

/** get the value string for attribute name
 *  if name is unfound in cur you get an empty string back
 *  this is the same behavior XMLAttrString provides. Without the complication
 *  of a composite type you don't need.
 */
std::string getXMLAttributeValue(const xmlNodePtr cur, const std::string_view name);

#endif
