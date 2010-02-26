//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef ESTOOLS_OHMMS_LIBXML_DEF_H
#define ESTOOLS_OHMMS_LIBXML_DEF_H
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <iosfwd>
#include <string>
#include <vector>
#include <algorithm>

template<typename _CharT>
inline void getNodeName(std::basic_string<_CharT>& cname, xmlNodePtr cur) {
  cname=(const char*)cur->name;
  for(int i=0; i<cname.size(); i++) cname[i]=tolower(cname[i]);
  //std::transform(cname.begin(), cname.end(), cname.begin(), std::tolower);
}

/**\file libxmldefs.h
 *\brief A collection of put/get functions to read from or write to a xmlNode defined in libxml2.
 *
 */
/*!\brief assign a value from a node. Use specialization for classes. 
 *\param a reference to a value to be assigned
 *\param cur current node
 *\return ture if successful
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
bool putContent(T& a, xmlNodePtr cur) {
  std::istringstream 
    stream((const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  return stream >> a;
}

template<class IT>
bool putContent(IT first, IT last, xmlNodePtr cur) {
  std::istringstream 
    stream((const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  bool success=true;
  while(success && first!=last) {
    success=(stream >> *first++);
  }
  return success;
}

/*!\fn bool getContent(const T& a, xmlNodePtr cur)
 *\brief write a value to a node. 
 *\param a reference to a value to be copied to a node
 *\param cur current node to which a content is copied
 *\return ture if successful
 *
 *Use specialization for classes. If operator << is implemented for a
 *class, no specialization is required. For instance, no
 *specialization is necessary for intrinsic data types.
  */
template<class T>
bool getContent(const T& a, xmlNodePtr cur) {
  std::stringstream s;
  s.setf(std::ios::scientific, std::ios::floatfield);
  s.precision(10);
  s << a;
  xmlNodeSetContent(cur,(const xmlChar*)(s.str().c_str()));
  return true;
}

/** assign vector<T> from a node. Create a temporary vector and make assignment.
 *\param a reference vector<T> 
 *\param cur current node to which a content is copied
 *\return ture if successful
 *
 *Specialization for std::vector<T> with data types T with operator >>,
 *e.g., std::vector<double>
 */
template<class T>
inline bool 
putContent(std::vector<T>& a, xmlNodePtr cur){
  std::istringstream 
    stream((const char*)
	   (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  std::vector<T> b;
  T t;
  while(!stream.eof()){ if(stream >> t) b.push_back(t);}
  a = b;
  return true;
}

/** write std::vector<T> to node. Each element is separated by a space.
 *\param a reference vector<T> 
 *\param cur current node to which a content is copied
 *\return ture if successful
 *
 *Specialization for std::vector<T> with data types T with operator <<.
 *This function is only for testing/debugging and will not perform
 *well for large vectors. Use HDF5 library for the real data.
 */
template<class T>
inline bool 
getContent(const std::vector<T>& a, xmlNodePtr cur) {
  std::stringstream s;
  s.precision(10);
  for(int i=0; i<a.size(); i++) s << ' ' <<  a[i];
  //xmlNodeAddContent(cur,(const xmlChar*)(s.str().c_str()));
  xmlNodeSetContent(cur,(const xmlChar*)(s.str().c_str()));
  return true;
}


#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
