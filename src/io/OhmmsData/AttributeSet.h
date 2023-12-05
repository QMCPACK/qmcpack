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
  std::map<std::string, std::unique_ptr<OhmmsElementBase>> m_param;

  bool get(std::ostream& os) const
  {
    for (const auto& [name, param] : m_param)
      param->get(os);
    return true;
  }

  /** add a new attribute
   *@param aparam reference the object which this attribute is assigned to.
   *@param aname the name of the added attribute
   *@param candidate_values candidate values to be checked against, the first element is the default value
   *@param status Tag status, See OhmmsParameter.h for more details
   */
  template<class PDT>
  void add(PDT& aparam,
           const std::string& aname,
           std::vector<PDT> candidate_values = {},
           TagStatus status                    = TagStatus::OPTIONAL)
  {
    if (auto it = m_param.find(aname); it == m_param.end())
      m_param[aname] = std::make_unique<OhmmsParameter<PDT>>(aparam, aname, std::move(candidate_values), status);
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
      if (auto it = m_param.find(aname); it != m_param.end())
      {
        std::istringstream stream((const char*)(att->children->content));
        it->second->put(stream);
      }
      att = att->next;
    }
    return true;
  }
};
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
