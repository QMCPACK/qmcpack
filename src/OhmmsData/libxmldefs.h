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
#include <typeinfo>
#include "OhmmsData/XMLParsingString.h"

template<typename _CharT>
inline void getNodeName(std::basic_string<_CharT>& cname, xmlNodePtr cur)
{
  cname = (const char*)cur->name;
  for (int i = 0; i < cname.size(); i++)
    cname[i] = tolower(cname[i]);
  //std::transform(cname.begin(), cname.end(), cname.begin(), std::tolower);
}

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
  std::istringstream stream(XMLNodeString{cur});
  std::vector<T> b;
  T t;
  while (!stream.eof())
  {
    if (stream >> t)
      b.push_back(t);
    else if (!stream.eof() && stream.fail())
    {
      std::cerr << "failbit detected when reading type (type_info::name) "
                << typeid(T).name() << ", value " << t << std::endl;
      stream.clear();
    }
  }
  a = b;
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


#endif
