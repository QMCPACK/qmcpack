//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file XMLString.h
 *
 * convert xmlNode contents into a std::string
 */
#ifndef QMCPLUSPLUS_XMLSTRING_H
#define QMCPLUSPLUS_XMLSTRING_H

#include <string>
#include <libxml/xmlmemory.h>
#include <libxml/tree.h>

class XMLString : public std::string
{
public:
  /// construct a string from an xmlNode
  XMLString(const xmlNodePtr cur)
  {
    xmlChar* node_char = xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1);
    if(node_char)
    {
      assign((const char*)node_char);
      xmlFree(node_char);
    }
  }

  /// expose base class constructors
  XMLString(const std::string& in) : std::string(in) { }

  XMLString(const char* in) : std::string(in) { }

  /// write a string to an xmlNode
  void setXMLNodeContent(xmlNodePtr cur) const
  {
    xmlNodeSetContent(cur, (const xmlChar*)(this->c_str()));
  }
};
#endif
