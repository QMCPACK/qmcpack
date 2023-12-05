//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Refactored from: ParameterSet.h
//////////////////////////////////////////////////////////////////////////////////////

#include "ParameterSet.h"
#include <vector>
#include <libxml/xmlmemory.h>
#include <libxml/tree.h>
#include "ModernStringUtils.hpp"
#include "libxmldefs.h"

bool ParameterSet::put(xmlNodePtr cur)
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

template<class PDT>
void ParameterSet::add(PDT& aparam, const std::string& aname_in, std::vector<PDT> candidate_values, TagStatus status)
{
  using namespace qmcplusplus;
  std::string aname(lowerCase(aname_in));
  if (auto it = m_param.find(aname); it == m_param.end())
  {
    m_param[aname] = std::make_unique<OhmmsParameter<PDT>>(aparam, aname, std::move(candidate_values), status);
  }
}

template<class PDT>
void ParameterSet::setValue(const std::string& aname_in, PDT aval)
{
  using namespace qmcplusplus;
  std::string aname(lowerCase(aname_in));
  if (auto it = m_param.find(aname); it != m_param.end())
  {
    (dynamic_cast<OhmmsParameter<PDT>&>(*it->second)).setValue(aval);
  }
}

template void ParameterSet::add(std::string&, const std::string&, std::vector<std::string>, TagStatus);
template void ParameterSet::add(qmcplusplus::astring&,
                                const std::string&,
                                std::vector<qmcplusplus::astring>,
                                TagStatus);
template void ParameterSet::add<int>(int&, const std::string&, std::vector<int>, TagStatus);
template void ParameterSet::add<bool>(bool&, const std::string&, std::vector<bool>, TagStatus);
template void ParameterSet::add<double>(double&, const std::string&, std::vector<double>, TagStatus);
template void ParameterSet::add<float>(float&, const std::string&, std::vector<float>, TagStatus);
template void ParameterSet::add<std::complex<double>>(std::complex<double>&,
                                                      const std::string&,
                                                      std::vector<std::complex<double>>,
                                                      TagStatus);
template void ParameterSet::add<std::complex<float>>(std::complex<float>&,
                                                     const std::string&,
                                                     std::vector<std::complex<float>>,
                                                     TagStatus);

template void ParameterSet::add<qmcplusplus::TinyVector<int, 3u>>(qmcplusplus::TinyVector<int, 3u>&,
                                                                  const std::string&,
                                                                  std::vector<qmcplusplus::TinyVector<int, 3u>>,
                                                                  TagStatus);

template void ParameterSet::setValue<int>(const std::string& aname_in, int aval);
