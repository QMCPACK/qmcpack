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

/** class to handle a set of parameters
 *
 *This may become an inherited class from OhmmsElementBase.
 */
struct ParameterSet : public OhmmsElementBase
{
  //  public std::map<std::string, OhmmsElementBase*> {

  typedef std::map<std::string, std::unique_ptr<OhmmsElementBase>> Container_t;
  typedef Container_t::iterator iterator;
  typedef Container_t::const_iterator const_iterator;

  Container_t m_param;

  ParameterSet(const char* aname = "parameter") : OhmmsElementBase(aname) {}

  inline bool get(std::ostream& os) const override
  {
    for (const auto& it : m_param)
      it.second->get(os);
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
    if (cur == NULL)
      return true; //handle empty node
    cur            = cur->xmlChildrenNode;
    bool something = false;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      tolower(cname);
      iterator it_tag = m_param.find(cname);
      if (it_tag == m_param.end())
      {
        if (cname == myName)
        {
          XMLAttrString aname(cur, "name");
          if (!aname.empty())
          {
            tolower(aname);
            iterator it = m_param.find(aname);
            if (it != m_param.end())
            {
              something = true;
              (*it).second->put(cur);
            }
          }
        }
      }
      else
      {
        something = true;
        (*it_tag).second->put(cur);
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
    std::string aname(aname_in);
    tolower(aname);
    iterator it = m_param.find(aname);
    if (it == m_param.end())
    {
      m_param[aname] = std::make_unique<OhmmsParameter<PDT>>(aparam, aname, std::move(candidate_values), status);
    }
  }

  template<class PDT>
  inline void setValue(const std::string& aname_in, PDT aval)
  {
    std::string aname(aname_in);
    tolower(aname);
    iterator it = m_param.find(aname);
    if (it != m_param.end())
    {
      (dynamic_cast<OhmmsParameter<PDT>*>((*it).second.get()))->setValue(aval);
    }
  }
};
#endif /*OHMMS_OHMMSPARAMETERSET_H*/
