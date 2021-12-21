//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_OHMMSPARAMETERSET_H
#define OHMMS_OHMMSPARAMETERSET_H

#include <map>
#include <string>
#include "OhmmsData/OhmmsParameter.h"
#include "ModernStringUtils.hpp"
/** class to handle a set of parameters
 *
 *This may become an inherited class from OhmmsElementBase.
 */
struct ParameterSet : public OhmmsElementBase
{
  //  public std::map<std::string, OhmmsElementBase*> {

  std::map<std::string, std::unique_ptr<OhmmsElementBase>> m_param;

  ParameterSet(const char* aname = "parameter") : OhmmsElementBase(aname) {}

  inline bool get(std::ostream& os) const override
  {
    for (const auto& [name, param] : m_param)
      param->get(os);
    return true;
  }

  inline bool put(std::istream& is) override { return true; }

  /** assign parameters to the set
   * @param cur the xml node to work on
   * @return true, if any valid parameter is processed.
   *
   * Accept both
   * - <aname> value </aname>
   * - <parameter name="aname"> value </parameter>
   * aname is converted into lower cases.
   */
  inline bool put(xmlNodePtr cur) override
  {
    using namespace qmcplusplus;
    if (cur == NULL)
      return true; //handle empty node
    cur            = cur->xmlChildrenNode;
    bool something = false;
    while (cur != NULL)
    {
      std::string cname(lowerCase(castXMLCharToChar(cur->name)));
      if (auto it_tag = m_param.find(cname); it_tag == m_param.end())
      {
        if (cname == myName)
        {
          std::string aname(lowerCase(getXMLAttributeValue(cur, "name")));
          if (!aname.empty())
          {
            if (auto it = m_param.find(aname); it != m_param.end())
            {
              something = true;
              it->second->put(cur);
            }
          }
        }
      }
      else
      {
        something = true;
        it_tag->second->put(cur);
      }
      cur = cur->next;
    }
    return something;
  }

  inline void reset() override {}

  /** add a new parameter corresponding to an xmlNode <parameter/>
   *@param aparam reference the object which this parameter is assigned to.
   *@param aname_in the value of the name attribute
   *@param candidate_values candidate values to be checked against, the first element is the default value
   *@param status Tag status, See OhmmsParameter.h for more details
   */
  template<class PDT>
  inline void add(PDT& aparam,
                  const std::string& aname_in,
                  std::vector<PDT>&& candidate_values = {},
                  TagStatus status                    = TagStatus::OPTIONAL)
  {
    using namespace qmcplusplus;
    std::string aname(lowerCase(aname_in));
    if (auto it = m_param.find(aname); it == m_param.end())
    {
      m_param[aname] = std::make_unique<OhmmsParameter<PDT>>(aparam, aname, std::move(candidate_values), status);
    }
  }

  template<class PDT>
  inline void setValue(const std::string& aname_in, PDT aval)
  {
    using namespace qmcplusplus;
    std::string aname(lowerCase(aname_in));
    if (auto it = m_param.find(aname); it != m_param.end())
    {
      (dynamic_cast<OhmmsParameter<PDT>&>(*it->second)).setValue(aval);
    }
  }
};
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
