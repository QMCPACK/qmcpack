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
#ifndef OHMMS_OHMMSPARAMETERSET_H
#define OHMMS_OHMMSPARAMETERSET_H

#include <map>
#include <string>
#include "OhmmsData/OhmmsParameter.h"

/** class to handle a set of parameters
 *
 *This may become an inherited class from OhmmsElementBase.
 */
struct ParameterSet {

  /** container for parameters that belong to this set
   */
  std::map<std::string,OhmmsElementBase*>  thisSet;

  bool get(std::ostream& os) {
    std::map<std::string,OhmmsElementBase*>::iterator it = thisSet.begin();
    while(it != thisSet.end()) {
      (*it).second->get(os);it++;
    }
    return true;
  }

  /** add a new parameter corresponding to an xmlNode <parameter/>
   *@param aparam reference the object which this parameter is assigned to.
   *@param aname the value of the name attribute
   *@param uname the value of the condition attribute
   *
   *The attributes of a parameter are
   * - name, the name of the parameter
   * - condition, the unit of the parameter
   *The condition will be used to convert the external unit to the internal unit.
   */
  template<class PDT>
  inline void add(PDT& aparam, const char* aname, const char* uname) {
    map<string,OhmmsElementBase*>::iterator it = thisSet.find(aname);
    if(it == thisSet.end()) {
      thisSet[aname] = new OhmmsParameter<PDT>(aparam,aname,uname);
    }
  }

  /** assign parameters to the set **/
  inline bool put(xmlNodePtr cur) {
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "parameter") {
	string aname = (const char*)(xmlGetProp(cur, (const xmlChar *) "name"));
	map<string,OhmmsElementBase*>::iterator it = thisSet.find(aname);
	if(it != thisSet.end()) {
	  (*it).second->put(cur);
	} 
      }
      cur=cur->next;
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
