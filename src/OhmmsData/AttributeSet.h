//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_XMLATTRIBUTESET_H
#define OHMMS_XMLATTRIBUTESET_H

#include <map>
#include <string>
#include "OhmmsData/OhmmsParameter.h"

/** class to handle a set of parameters
 *
 *This may become an inherited class from OhmmsElementBase.
 */
struct OhmmsAttributeSet
{
  typedef std::map<std::string, OhmmsElementBase*>  Container_t;
  typedef Container_t::iterator iterator;
  typedef Container_t::const_iterator const_iterator;

  xmlNodePtr myNode;
  Container_t m_param;

  inline OhmmsAttributeSet(): myNode(0) {}

  ~OhmmsAttributeSet()
  {
    iterator it(m_param.begin());
    iterator it_end(m_param.end());
    while(it!=it_end)
    {
      delete (*it).second;
      ++it;
    }
  }

  bool get(std::ostream& os) const
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
  void add(PDT& aparam, const std::string& aname)
  {
    iterator it(m_param.find(aname));
    if(it == m_param.end())
    {
      m_param[aname] = new OhmmsParameter<PDT>(aparam,aname.c_str(),"none");
    }
  }

  /** assign parameters to the set
   *@param cur the xml node to work on
   *@return true, if any valid parameter is processed.
   */
  bool put(xmlNodePtr cur)
  {
    xmlAttrPtr att = cur->properties;
    while(att != NULL)
    {
      std::string aname((const char*)(att->name));
      iterator it = m_param.find(aname);
      if(it != m_param.end())
      {
        std::istringstream stream((const char*)(att->children->content));
        (*it).second->put( stream);
      }
      att=att->next;
    }
    return true;
  }
};
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
