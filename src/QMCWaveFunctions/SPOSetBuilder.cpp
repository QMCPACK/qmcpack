//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <QMCWaveFunctions/SPOSetBuilder.h>
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
SPOSetBuilder::SPOSetBuilder(const std::string& SPO_type_name_in, Communicate* comm)
    : MPIObjectBase(comm), legacy(true), SPO_type_name(SPO_type_name_in)
{
  reserve_states();
}


void SPOSetBuilder::reserve_states(int nsets)
{
  int sets_needed = nsets - states.size();
  if (sets_needed > 0)
    for (int s = 0; s < sets_needed; ++s)
      states.push_back(new SPOSetInfo());
}


SPOSet* SPOSetBuilder::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info)
{
  APP_ABORT("BasisSetBase::createSPOSet(cur,input_info) has not been implemented");
  return 0;
}


SPOSet* SPOSetBuilder::createSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name("");
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.put(cur);

  if (spo_object_name.empty())
    app_warning() << "SPOSet object name not given in the input!" << std::endl;

  app_summary() << std::endl;
  app_summary() << "   Single particle orbitals (SPO)" << std::endl;
  app_summary() << "   ------------------------------" << std::endl;
  app_summary() << "    Name: " << spo_object_name << "   Type: " << SPO_type_name
                << "   Builder class name: " << ClassName << std::endl;
  app_summary() << std::endl;

  // read specialized sposet construction requests
  //   and translate them into a set of orbital indices
  SPOSetInputInfo input_info(cur);

  // process general sposet construction requests
  //   and preserve legacy interface
  SPOSet* sposet = 0;
  if (legacy && input_info.legacy_request)
    sposet = createSPOSetFromXML(cur);
  else
    sposet = createSPOSet(cur, input_info);

  // remember created sposets
  if (sposet)
  {
    //sposet->put(cur); //initialize C and other internal containers
    sposets.push_back(sposet);
  }
  else
    APP_ABORT("SPOSetBuilder::createSPOSet  sposet creation failed");

  if (!spo_object_name.empty() && sposet->getName().empty())
    sposet->setName(spo_object_name);
  if (sposet->getName().empty())
    app_warning() << "SPOSet object doesn't have a name." << std::endl;
  if (!spo_object_name.empty() && sposet->getName() != spo_object_name)
    app_warning() << "SPOSet object name mismatched! input name: " << spo_object_name
                  << "   object name: " << sposet->getName() << std::endl;

  return sposet;
}

} // namespace qmcplusplus
