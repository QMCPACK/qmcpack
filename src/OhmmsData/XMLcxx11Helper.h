//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 QMCPACK developers.
//
// File developed by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California Berkeley
//
// File created by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef OHMMS_XMLCXX11HELPER_H
#define OHMMS_XMLCXX11HELPER_H

#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  gets the supplied xml node's children in std::list<xmlNodePtr> form, to facilitate
///         convenient looping like:
///
///           for ( auto ptr : getXMLChildList(cur) ) {
///             // Do stuff with ptr, which iterates over the children
///           }
///
/// \param[in,out]  cur      pointer to the xml node whose children we are interested in
///
///////////////////////////////////////////////////////////////////////////////////////////////////
inline std::list<xmlNodePtr> getXMLChildList(xmlNodePtr cur) {

  // initialize list to return
  std::list<xmlNodePtr> retval;

  // fill list with the supplied node's children
  for (xmlNodePtr chl = cur->xmlChildrenNode; chl != NULL; chl = chl->next)
    retval.push_back(chl);

  // return the list
  return retval;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns a list in std::list<xmlNodePtr> form of all the supplied xml node's children
///         that have the requested name.
///         Facilitates convenient looping like this:
///
///           for ( auto ptr : getXMLChildList(cur, "childName") ) {
///             // Do stuff with ptr, which iterates over the children
///             // that have name "childName".
///           }
///
/// \param[in,out]  cur      pointer to the xml node whose children we are interested in
/// \param[in]      name     the name we are looking for
///
///////////////////////////////////////////////////////////////////////////////////////////////////
inline std::list<xmlNodePtr> getXMLChildList(xmlNodePtr cur, const std::string & name) {

  // initialize list to return
  std::list<xmlNodePtr> retval;

  // loop over the children, adding those with the correct name
  for ( auto chl : getXMLChildList(cur) )
    if ( std::string((const char*)chl->name) == name )
      retval.push_back(chl);

  // return the list
  return retval;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Termination step of the recursion when no arguments are left.
///
/// \param[in,out]  att      attribute set object
///
///////////////////////////////////////////////////////////////////////////////////////////////////
inline void getXMLAttributesRecurse(OhmmsAttributeSet & att) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Propagation step of the recursion that prepares to read the next attribute value
///
/// \param[in,out]  att      attribute set object
/// \param[out]     arg      the variable to store the attribute value in
/// \param[in]      name     the name of the attribute to read the value from
/// \param[in,out]  args     variable-length argument list of remaining variables and their names
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, typename... Args>
inline void getXMLAttributesRecurse(OhmmsAttributeSet & att, T & arg, const std::string & name, Args & ... args) {

  // set up to read the next attribute value
  att.add(arg, name);

  // process any remaining arguments (or terminate recursion if there are none left)
  getXMLAttributesRecurse(att, args...);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Reads attributes of the supplied xml node and records those named in the
///         variable-length argument list in the supplied variables.
///         Use like this:
///
///           // cur is a pointer to an xml node
///
///           double x;
///           int a;
///           std::string s;
///           getXMLAttributes(cur, x, "xname", a, "aname", s, "sname");
///
///           // Now the value of the attribute named xname is stored in x,
///           // the value of the attribute named aname is stored in a,
///           // and the value of the attribute named sname is stored in s.
///
/// \param[in,out]  p        pointer to the xml node whose attributes we are interested in
/// \param[in,out]  args     variable-length argument list of variables and their names
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<typename... Args>
inline bool getXMLAttributes(xmlNodePtr & p, Args & ... args) {

  // check if argument list is a sane length
  if ( ( sizeof...(args) ) % 2 != 0 )
    APP_ABORT("getXMLAttributes requires an even number of template arguments");

  // read in attribute values if any have been requested
  if ( ( sizeof...(args) ) > 0 ) {

    // prepare attribute set object used for reading
    OhmmsAttributeSet att;

    // recursively process arguments to find out which attributes to read
    getXMLAttributesRecurse(att, args...);

    // read the attributes into the variables
    return att.put(p);

  }

  // return true if we were given nothing to read
  return true;

}

}

#endif
