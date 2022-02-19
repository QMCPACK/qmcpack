//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file ExampleHeBuilder.cpp
 *@brief Example builder for simple He wavefunction
 */
#include "ExampleHeBuilder.h"
#include "QMCWaveFunctions/ExampleHeComponent.h"
#include "OhmmsData/AttributeSet.h"


namespace qmcplusplus
{
ExampleHeBuilder::ExampleHeBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets)
    : WaveFunctionComponentBuilder(comm, p), ptclPool(psets), els(p)
{}

std::unique_ptr<WaveFunctionComponent> ExampleHeBuilder::buildComponent(xmlNodePtr cur)
{
  std::string ion_name = "ion0";
  OhmmsAttributeSet oAttrib;
  oAttrib.add(ion_name, "source");
  oAttrib.put(cur);

  auto ion_it = ptclPool.find(ion_name);
  if (ion_it == ptclPool.end())
  {
    app_error() << " Ion particle set not found  = " << ion_name << std::endl;
    APP_ABORT("Ion not found");
  }
  auto WF = std::make_unique<ExampleHeComponent>(*(ion_it->second), els);
  WF->put(cur);
  return WF;
}

} // namespace qmcplusplus
