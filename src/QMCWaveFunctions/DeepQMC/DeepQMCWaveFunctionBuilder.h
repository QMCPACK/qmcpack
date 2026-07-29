//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DEEPQMCWAVEFUNCTIONBUILDER_H
#define QMCPLUSPLUS_DEEPQMCWAVEFUNCTIONBUILDER_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{

class DeepQMCWaveFunctionBuilder : public WaveFunctionComponentBuilder
{
public:
  DeepQMCWaveFunctionBuilder(Communicate* comm, ParticleSet& target, const PSetMap& psets);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  const PSetMap& ptcl_pool_;
};

} // namespace qmcplusplus

#endif
