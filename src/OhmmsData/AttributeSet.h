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

/** class to handle a set of attributes of an xmlNode
 */
struct OhmmsAttributeSet
{
  std::map<std::string, OhmmsElementBase*> m_param;

  ~OhmmsAttributeSet()
  {
    
    decltype(m_param)::iterator it(m_param.begin());
    decltype(m_param)::iterator it_end(m_param.end());
    while (it != it_end)
    {
      delete (*it).second;
      ++it;
    }
  }

  bool get(std::ostream& os) const
  {
    decltype(m_param)::const_iterator it(m_param.begin());
    decltype(m_param)::const_iterator it_end(m_param.end());
    while (it != it_end)
    {
      (*it).second->get(os);
      ++it;
    }
    return true;
  }

  /** add a new attribute
   *@param aparam reference the object which this attribute is assigned to.
   *@param aname the name of the added attribute
   */
  template<class PDT>
  void add(PDT& aparam, const std::string& aname)
  {
    decltype(m_param)::iterator it(m_param.find(aname));
    if (it == m_param.end())
    {
      m_param[aname] = new OhmmsParameter<PDT>(aparam, aname.c_str(), "none");
    }
  }

  /** assign attributes to the set
   *@param cur the xml node to work on
   *@return true, if any valid parameter is processed.
   */
  bool put(xmlNodePtr cur)
  {
    xmlAttrPtr att = cur->properties;
    while (att != NULL)
    {
      std::string aname((const char*)(att->name));
      decltype(m_param)::iterator it = m_param.find(aname);
      if (it != m_param.end())
      {
        std::istringstream stream((const char*)(att->children->content));
        (*it).second->put(stream);
      }
      att = att->next;
    }
    return true;
  }
};
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
