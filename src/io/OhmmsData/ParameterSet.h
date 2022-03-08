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
#include <complex>
#include "OhmmsData/OhmmsParameter.h"
#include "ModernStringUtils.hpp"
#include "OhmmsPETE/TinyVector.h"

/** class to handle a set of parameters
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
  bool put(xmlNodePtr cur) override;

  inline void reset() override {}

  /** add a new parameter corresponding to an xmlNode <parameter/>
   *@param aparam reference the object which this parameter is assigned to.
   *@param aname_in the value of the name attribute
   *@param candidate_values candidate values to be checked against, the first element is the default value
   *@param status Tag status, See OhmmsParameter.h for more details
   */
  template<class PDT>
  void add(PDT& aparam,
           const std::string& aname_in,
           std::vector<PDT> candidate_values = {},
           TagStatus status                    = TagStatus::OPTIONAL);

  template<class PDT>
  void setValue(const std::string& aname_in, PDT aval);
};

extern template void ParameterSet::add<std::string>(std::string&,
                                                    const std::string&,
                                                    std::vector<std::string>,
                                                    TagStatus);
extern template void ParameterSet::add<bool>(bool&, const std::string&, std::vector<bool>, TagStatus);
extern template void ParameterSet::add<int>(int&, const std::string&, std::vector<int>, TagStatus);
extern template void ParameterSet::add<double>(double&, const std::string&, std::vector<double>, TagStatus);
extern template void ParameterSet::add<float>(float&, const std::string&, std::vector<float>, TagStatus);
extern template void ParameterSet::add<std::complex<double>>(std::complex<double>&,
                                                             const std::string&,
                                                             std::vector<std::complex<double>>,
                                                             TagStatus);
extern template void ParameterSet::add<std::complex<float>>(std::complex<float>&,
                                                            const std::string&,
                                                            std::vector<std::complex<float>>,
                                                            TagStatus);
extern template void ParameterSet::add<qmcplusplus::TinyVector<int, 3u>>(
    qmcplusplus::TinyVector<int, 3u>&,
    const std::string&,
    std::vector<qmcplusplus::TinyVector<int, 3u>>,
    TagStatus);

extern template void ParameterSet::setValue<int>(const std::string& aname_in, int aval);

#endif /*OHMMS_OHMMSPARAMETERSET_H*/
