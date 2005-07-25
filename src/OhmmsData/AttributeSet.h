//////////////////////////////////////////////////////////////////
// (c) Copyright 2004-  by Jeongnim Kim
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
#ifndef OHMMS_XMLATTRIBUTESET_H
#define OHMMS_XMLATTRIBUTESET_H

#include <map>
#include <string>
#include "OhmmsData/OhmmsParameter.h"

/** class to handle a set of parameters
 *
 *This may become an inherited class from OhmmsElementBase.
 */
struct OhmmsAttributeSet: public std::map<std::string, OhmmsElementBase*> 
{

  xmlNodePtr myNode;

  inline OhmmsAttributeSet(): myNode(0) {}

  ~OhmmsAttributeSet() {
    iterator it = begin();
    iterator it_end = end();
    while(it!=it_end) {delete (*it).second; ++it;}
  }

  bool get(std::ostream& os) const {
    const_iterator it = begin();
    const_iterator it_end = end();
    while(it != it_end) {
      (*it).second->get(os);++it;
    }
    return true;
  }

  /** add a new parameter corresponding to an xmlNode <parameter/>
   *@param aparam reference the object which this parameter is assigned to.
   *@param aname the value of the name attribute
   *
   *The attributes of a parameter are
   * - name, the name of the parameter
   * - condition, the unit of the parameter
   *The condition will be used to convert the external unit to the internal unit.
   */
  template<class PDT>
  inline void add(PDT& aparam, const string& aname) {
    iterator it = find(aname);
    if(it == end()) {
      operator[](aname) = new OhmmsParameter<PDT>(aparam,aname.c_str(),"none");
    }
  }

  /** assign parameters to the set 
   *@param cur the xml node to work on
   *@return true, if any valid parameter is processed.
   */
  inline bool put(xmlNodePtr cur) {
    xmlAttrPtr att = cur->properties;
    while(att != NULL) {
      string aname((const char*)(att->name));
      iterator it = find(aname);
      if(it != end()) {
        std::istringstream stream((const char*)(att->children->content));
        (*it).second->put(stream);
      } 
      att=att->next;
    }
    return true;
  }
}; 
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
