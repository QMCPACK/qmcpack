//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//
// File created by: Brett Van Der Goetz, bvdg@berkeley.edu, University of California at Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COUNTING_JASTROW_BUILDER_H
#define QMCPLUSPLUS_COUNTING_JASTROW_BUILDER_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
class CountingJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  // voronoi constructor
  CountingJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source);
  // normalized gaussian constructor
  CountingJastrowBuilder(Communicate* comm, ParticleSet& target);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  std::string NameOpt;
  std::string TypeOpt;
  std::string RegionOpt;
  std::string SourceOpt;

  ParticleSet* SourcePtcl;

  std::unique_ptr<WaveFunctionComponent> createCJ(xmlNodePtr cur);
};

} // namespace qmcplusplus

#endif
