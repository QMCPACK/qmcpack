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


/** @file XMLParsingString.h
 *
 * XMLNodeString convert xmlNode contents into a std::string
 * XMLAttrString convert one xmlNode attribute into a std::string
 */
#ifndef QMCPLUSPLUS_XMLSTRING_H
#define QMCPLUSPLUS_XMLSTRING_H

#include <string>
#include <libxml/xmlmemory.h>
#include <libxml/tree.h>

/** convert xmlNode contents into a std::string
 */
class XMLNodeString : public std::string
{
public:
  /// construct a string from an xmlNode
  XMLNodeString(const xmlNodePtr cur)
  {
    xmlChar* node_char = xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1);
    if(node_char)
    {
      assign((const char*)node_char);
      xmlFree(node_char);
    }
  }

  /// expose base class constructors
  XMLNodeString(const std::string& in) : std::string(in) { }

  XMLNodeString(const char* in) : std::string(in) { }

  /// write a string to an xmlNode
  void setXMLNodeContent(xmlNodePtr cur) const
  {
    xmlNodeSetContent(cur, (const xmlChar*)(this->c_str()));
  }
};

/** convert an xmlNode attribute into a std::string
 * if parsing multiple attributes is needed, use OhmmsAttributeSet
 */
class XMLAttrString : public std::string
{
  const std::string attribute_name_;

public:
  /// construct a string from an xmlNode
  XMLAttrString(const xmlNodePtr cur, const char* name)
      : attribute_name_(name)
  {
    xmlChar* attr_char = xmlGetProp(cur, (const xmlChar*)name);
    if(attr_char)
    {
      assign((const char*)attr_char);
      xmlFree(attr_char);
    }
  }

  /// expose base class constructors
  XMLAttrString(const std::string& in, const std::string& name) : std::string(in), attribute_name_(name) { }

  XMLAttrString(const char* in, const char* name) : std::string(in), attribute_name_(name) { }

  /// write a string to an xmlNode
  void createXMLAttribute(xmlNodePtr cur) const
  {
    xmlNewProp(cur, (const xmlChar*)attribute_name_.c_str(), (const xmlChar*)this->c_str());
  }

  /// write a string to an xmlNode
  void setXMLAttribute(xmlNodePtr cur) const
  {
    xmlSetProp(cur, (const xmlChar*)attribute_name_.c_str(), (const xmlChar*)this->c_str());
  }
};
#endif
