//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef ESTOOLS_OHMMS_LIBXML_DEF_H
#define ESTOOLS_OHMMS_LIBXML_DEF_H
#include <iostream>
#include <sstream>
#include <iomanip>
#include <iosfwd>
#include <string>
#include <vector>
#include <algorithm>
#include "XMLParsingString.h"
#include "string_utils.h"

/**\file libxmldefs.h
 *\brief A collection of put/get functions to read from or write to a xmlNode defined in libxml2.
 *
 */
/*!\brief assign a value from a node. Use specialization for classes.
 *\param a reference to a value to be assigned
 *\param cur current node
 *\return true if successful
 *
 *If operator >> is implemented for a class, no specialization is required.
 *For instance, no specialization is necessary for intrinsic data types.
 *A simle class that reads the temperature from a xml node
 *
 \<parameter name="temperature" condition="K"\>100\</parameter\>
 *
 \code
  struct A: public OhmmsElementBase {

    double Temperature;

    bool put(xmlNodePtr cur) {
      putContent(Temperature,cur);
    }
  };
 \endcode
 */

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

std::string getNodeName(xmlNodePtr cur);

/** replaces a's value with the first "element" in the "string"
 *  returned by XMLNodeString{cur}.
 *  but if that string is empty then then doesn't touch a.
 *  See documentation of the interaction between operator>>
 *  streams and the FormattedInputFunction requirement in the c++ standard.
 */
template<class T>
bool putContent(T& a, xmlNodePtr cur)
{
  std::istringstream stream(XMLNodeString{cur});
  stream >> a;
  return !stream.fail();
}

template<class IT>
bool putContent(IT first, IT last, xmlNodePtr cur)
{
  std::istringstream stream(XMLNodeString{cur});
  bool success = true;
  while (success && first != last)
  {
    stream >> *first++;
    success = !stream.fail();
  }
  return success;
}

/*!\fn bool getContent(const T& a, xmlNodePtr cur)
 *\brief write a value to a node.
 *\param a reference to a value to be copied to a node
 *\param cur current node to which a content is copied
 *\return true if successful
 *
 *Use specialization for classes. If operator << is implemented for a
 *class, no specialization is required. For instance, no
 *specialization is necessary for intrinsic data types.
  */
template<class T>
bool getContent(const T& a, xmlNodePtr cur)
{
  if (cur->children == NULL)
    return false;
  std::stringstream s;
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(10);
  s << a;
  const XMLNodeString node_string(s.str());
  node_string.setXMLNodeContent(cur);
  return true;
}

/** assign std::vector<T> from a node. Create a temporary vector and make assignment.
 *\param a reference std::vector<T>
 *\param cur current node to which a content is copied
 *\return true if successful
 *
 *Specialization for std::vector<T> with data types T with operator >>,
 *e.g., std::vector<double>
 */
template<class T>
inline bool putContent(std::vector<T>& a, const xmlNodePtr cur)
{
  if (cur->children == NULL)
    return false;
  a = qmcplusplus::convertStrToVec<T>(XMLNodeString{cur});
  return true;
}

/** write std::vector<T> to node. Each element is separated by a space.
 *\param a reference std::vector<T>
 *\param cur current node to which a content is copied
 *\return true if successful
 *
 *Specialization for std::vector<T> with data types T with operator <<.
 *This function is only for testing/debugging and will not perform
 *well for large vectors. Use HDF5 library for the real data.
 */
template<class T>
inline bool getContent(const std::vector<T>& a, xmlNodePtr cur)
{
  std::stringstream s;
  s.precision(10);
  for (int i = 0; i < a.size(); i++)
    s << ' ' << a[i];
  const XMLNodeString node_string(s.str());
  node_string.setXMLNodeContent(cur);
  return true;
}

/** process through all the children of an XML element
 * F is a lambda or functor
 * void F/[](const std::string& cname, const xmlNodePtr element) {...}
 */
template<class F>
void processChildren(const xmlNodePtr cur, const F& functor)
{
  xmlNodePtr element = cur->children;
  while (element != NULL)
  {
    std::string cname((const char*)(element->name));
    functor(cname, element);
    element = element->next;
  }
}

#endif
