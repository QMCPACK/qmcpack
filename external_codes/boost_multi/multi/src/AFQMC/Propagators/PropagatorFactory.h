//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_PROPAGATORFACTORY_H
#define QMCPLUSPLUS_AFQMC_PROPAGATORFACTORY_H

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/RandomGenerator.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Propagators/Propagator.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class PropagatorFactory
{
public:
  PropagatorFactory(std::map<std::string, AFQMCInfo>& info) : InfoMap(info) {}

  ~PropagatorFactory() {}

  bool is_constructed(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in WavefunctionFactory::is_constructed(string&): Missing xml block. \n");
    auto p0 = propagators.find(ID);
    if (p0 == propagators.end())
      return false;
    else
      return true;
  }

  // returns a pointer to the base Propagator class associated with a given ID
  Propagator& getPropagator(TaskGroup_& TG, const std::string& ID, Wavefunction& wfn, RandomBase<RealType>& rng)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      APP_ABORT(" Error in PropagatorFactory::getPropagator(string&): Missing xml block. \n");
    auto p0 = propagators.find(ID);
    if (p0 == propagators.end())
    {
      auto newp = propagators.insert(std::make_pair(ID, buildPropagator(TG, xml->second, wfn, rng)));
      if (not newp.second)
        APP_ABORT(" Error: Problems building new propagator in PropagatorFactory::getPropagator(string&). \n");
      return (newp.first)->second;
    }
    else
      return p0->second;
  }

  xmlNodePtr getXML(const std::string& ID)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml == xmlBlocks.end())
      return nullptr;
    else
      return xml->second;
  }

  // adds a xml block from which a Propagator can be built
  void push(const std::string& ID, xmlNodePtr cur)
  {
    auto xml = xmlBlocks.find(ID);
    if (xml != xmlBlocks.end())
      APP_ABORT("Error: Repeated Propagator block in PropagatorFactory. Propagator names must be unique. \n");
    xmlBlocks.insert(std::make_pair(ID, cur));
  }

protected:
  // reference to container of AFQMCInfo objects
  std::map<std::string, AFQMCInfo>& InfoMap;

  // generates a new Propagator and returns the pointer to the base class
  Propagator buildPropagator(TaskGroup_& TG, xmlNodePtr cur, Wavefunction& wfn, RandomBase<RealType>& rng)
  {
    std::string type("afqmc");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(type, "type");
    oAttrib.put(cur);

    app_log() << "\n****************************************************\n"
              << "               Initializing Propagator \n"
              << "****************************************************\n"
              << std::endl;

    if (type == "afqmc")
      return buildAFQMCPropagator(TG, cur, wfn, rng);
    else
    {
      app_error() << "Unknown Propagator type in PropagatorFactory::buildPropagator(): " << type << std::endl;
      APP_ABORT(" Error: Unknown Propagator type in PropagatorFactory::buildPropagator() \n");
    }
    return Propagator{};
  }

  Propagator buildAFQMCPropagator(TaskGroup_& TG, xmlNodePtr cur, Wavefunction& wfn, RandomBase<RealType>& r);

  std::map<std::string, xmlNodePtr> xmlBlocks;

  std::map<std::string, Propagator> propagators;
};


} // namespace afqmc

} // namespace qmcplusplus

#endif
