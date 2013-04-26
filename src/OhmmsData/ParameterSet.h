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
struct ParameterSet: public OhmmsElementBase
{
//  public std::map<std::string, OhmmsElementBase*> {

  typedef std::map<std::string, OhmmsElementBase*>  Container_t;
  typedef Container_t::iterator       iterator;
  typedef Container_t::const_iterator const_iterator;

  Container_t m_param;

  ParameterSet(const char* aname="parameter"):
    OhmmsElementBase(aname)
  {
  }

  ~ParameterSet()
  {
    iterator it(m_param.begin());
    iterator it_end(m_param.end());
    while(it!=it_end)
    {
      delete (*it).second;
      ++it;
    }
  }

  inline bool get(std::ostream& os) const
  {
    const_iterator it(m_param.begin());
    const_iterator it_end(m_param.end());
    while(it != it_end)
    {
      (*it).second->get(os);
      ++it;
    }
    return true;
  }

  inline bool put(std::istream& is)
  {
    return true;
  }

  /** assign parameters to the set
   *@param cur the xml node to work on
   *@return true, if any valid parameter is processed.
   */
  inline bool put(xmlNodePtr cur)
  {
    if(cur == NULL)
      return true;//handle empty node
    cur = cur->xmlChildrenNode;
    bool something = false;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      iterator it_tag = m_param.find(cname);
      if(it_tag == m_param.end())
      {
        if(cname == myName)
        {
          const xmlChar* aptr= xmlGetProp(cur, (const xmlChar *) "name");
          if(aptr)
          {
            //string aname = (const char*)(xmlGetProp(cur, (const xmlChar *) "name"));
            iterator it = m_param.find((const char*)aptr);
            if(it != m_param.end())
            {
              something =true;
              (*it).second->put(cur);
            }
          }
        }
      }
      else
      {
        something =true;
        (*it_tag).second->put(cur);
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
// INLINE_ALL void add(PDT& aparam, const char* aname, const char* uname) {
  inline void add(PDT& aparam, const char* aname, const char* uname)
  {
    iterator it = m_param.find(aname);
    if(it == m_param.end())
    {
      m_param[aname] = new OhmmsParameter<PDT>(aparam,aname,uname);
    }
  }

  template<class PDT>
  //INLINE_ALL void setValue(const std::string& aname, PDT aval){
  inline void setValue(const std::string& aname, PDT aval)
  {
    iterator it = m_param.find(aname);
    if(it != m_param.end())
    {
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
