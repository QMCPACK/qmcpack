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
#include <string_view>
#include <libxml/xmlmemory.h>
#include <libxml/tree.h>

// I think these should be in libxmldefs.h
// but it includes XMLParsingString so these cast functions need to be here.

/** xmlChar* to char* cast
 *  certainly not typesafe but at this interface between libxml and c++
 *  well it beats c style casts and might be more correct than a reinterpret.
 *
 *  This is fine with UTF-8 bytes going into a std::string.
 */
inline char* castXMLCharToChar(xmlChar* c) { return static_cast<char*>(static_cast<void*>(c)); }

/** unsafe const xmlChar* to const char* cast
 */
inline const char* castXMLCharToChar(const xmlChar* c) { return static_cast<const char*>(static_cast<const void*>(c)); }

/** unsafe char* to xmlChar* cast
 */
inline xmlChar* castCharToXMLChar(char* c) { return static_cast<xmlChar*>(static_cast<void*>(c)); }

/** unsafe const char* to const xmlChar* cast
 */
inline const xmlChar* castCharToXMLChar(const char* c)
{
  return static_cast<const xmlChar*>(static_cast<const void*>(c));
}

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

/** convert an xmlNode attribute into a name value pair.
 * attribute_ is an empty string if attribute_name_
 is not found.
 * if parsing multiple attributes is needed, use OhmmsAttributeSet
 */
class XMLAttrString
{
  // ordering same as crazy value, name signature of constructor
  std::string attribute_;
  const std::string attribute_name_;

public:
  /// construct a string from an xmlNode
  XMLAttrString(const xmlNodePtr cur, const char* name) : attribute_name_(name)
  {
    xmlChar* attr_char = xmlGetProp(cur, castCharToXMLChar(name));
    if (attr_char)
    {
      attribute_.assign(castXMLCharToChar(attr_char));
      xmlFree(attr_char);
    }
  }

  /// expose base class constructors
  XMLAttrString(const std::string& in, const std::string& name) : attribute_(in), attribute_name_(name) {}

  XMLAttrString(const char* in, const char* name) : attribute_(in), attribute_name_(name) {}

  /** i.e. does it have a value
   */
  auto hasValue() const { return !attribute_.empty(); };

  /// write a string to an xmlNode
  void createXMLAttribute(xmlNodePtr cur) const
  {
    xmlNewProp(cur, castCharToXMLChar(attribute_name_.c_str()), castCharToXMLChar(attribute_.c_str()));
  }

  /// write a string to an xmlNode
  void setXMLAttribute(xmlNodePtr cur) const
  {
    xmlSetProp(cur, castCharToXMLChar(attribute_name_.c_str()), castCharToXMLChar(attribute_.c_str()));
  }

  /** retain surprising behavior where this inherently name value object 
   *  can also operate as just a string of its value
   *  \todo remove this its a footgun and makes for less readable code
   */
  //operator const std::string&() { return attribute_; }

  const std::string& getName() const { return attribute_name_; }
  const std::string& getValue() const { return attribute_; }
  void setValue(const std::string_view value) { attribute_.assign(value); }
};
#endif
