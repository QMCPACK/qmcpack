//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERSETFACTORY_H
#define QMCPLUSPLUS_AFQMC_WALKERSETFACTORY_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "OhmmsData/libxmldefs.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Walkers/WalkerSet.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class WalkerSetFactory
{
public:
  WalkerSetFactory(std::map<std::string, AFQMCInfo>& info) : InfoMap(info) {}

  ~WalkerSetFactory() {}

  bool is_constructed(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in WalkerSetFactory::is_constructed(string&): Missing xml block. \n");
    auto wlk = handlers.find(ID);
    if (wlk == handlers.end())
      return false;
    else
      return true;
  }

  WalkerSet& getWalkerSet(TaskGroup_& TG, const std::string& ID, RandomBase<RealType>& rng)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT("Error: Missing xml Block in WalkerSetFactory::getWalkerSet(string&). \n");
    auto wlk = handlers.find(ID);
    if (wlk == handlers.end())
    {
      auto newwlk = handlers.insert(std::make_pair(ID, buildHandler(TG, xml->second, rng)));
      if (!newwlk.second)
        APP_ABORT(" Error: Problems inserting new hamiltonian in WalkerSetFactory::getHandler(streing&). \n");
      return (newwlk.first)->second;
    }
    else
      return wlk->second;
  }

  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

  // adds a xml block from which a WalkerSet can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated WalkerSet block in WalkerSetFactory. WalkerSet names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID, cur));
  }

protected:
  // reference to container of AFQMCInfo objects
  std::map<std::string, AFQMCInfo>& InfoMap;

  // generates a new WalkerSet and returns the pointer to the base class
  WalkerSet buildHandler(TaskGroup_& TG, xmlNodePtr cur, RandomBase<RealType>& rng)
  {
    std::string type("shared");
    std::string info("info0");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(info, "info");
    oAttrib.add(type, "type");
    oAttrib.put(cur);

    if (InfoMap.find(info) == InfoMap.end())
    {
      app_error() << "ERROR: Undefined info:" << info << "  \n";
      APP_ABORT("");
    }

    // keep like this until you have another choice and a variant framework in place
    if (type != "shared")
      APP_ABORT(" Error: Unknown WalkerSet type in WalkerSetFactory::buildHandler(). \n");

    return WalkerSet(TG, cur, InfoMap[info], rng);
  }

  std::map<std::string, xmlNodePtr> xmlBlocks;

  std::map<std::string, WalkerSet> handlers;
};
} // namespace afqmc
} // namespace qmcplusplus

#endif
