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
struct ParameterSet: public OhmmsElementBase,
  public std::map<std::string, OhmmsElementBase*> {

  ParameterSet(const char* aname="parameter"):
    OhmmsElementBase(aname) {
  }

  ~ParameterSet() {
    iterator it = begin();
    iterator it_end = end();
    while(it!=it_end) { delete (*it).second; ++it;}
  }

  inline bool get(std::ostream& os) const {
    const_iterator it = begin();
    const_iterator it_end = end();
    while(it != it_end) {
      (*it).second->get(os);++it;
    }
    return true;
  }

  inline bool put(std::istream& is) {
    return true;
  }

  /** assign parameters to the set 
   *@param cur the xml node to work on
   *@return true, if any valid parameter is processed.
   */
  inline bool put(xmlNodePtr cur) {
    cur = cur->xmlChildrenNode;
    bool something = false;
    while(cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname == myName) {
        const xmlChar* aptr= xmlGetProp(cur, (const xmlChar *) "name");
        if(aptr) {
	  //string aname = (const char*)(xmlGetProp(cur, (const xmlChar *) "name"));
	  iterator it = find((const char*)aptr);
	  if(it != end()) {
            something =true;
	    (*it).second->put(cur);
	  } 
        }
      }
      cur=cur->next;
    }
    return something;
  }

  inline void reset() { }

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
    iterator it = find(aname);
    if(it == end()) {
      operator[](aname) = new OhmmsParameter<PDT>(aparam,aname,uname);
    }
  }

  template<class PDT>
  inline void setValue(const std::string& aname, PDT aval){
    iterator it = find(aname);
    if(it != end()) {
       (dynamic_cast<OhmmsParameter<PDT>*>((*it).second))->setValue(aval);
    }
  }

}; 
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
