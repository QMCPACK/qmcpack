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
struct CountingJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  // voronoi constructor
  CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet& source);
  // normalized gaussian constructor
  CountingJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);

private:
  ///jastrow/@name
  std::string NameOpt;
  ///jastrow/@type
  std::string TypeOpt;
  ///jastrow/@region
  std::string RegionOpt;
  ///jastrow/@source
  std::string SourceOpt;

  ParticleSet* SourcePtcl;

  bool createCJ(xmlNodePtr cur);
};

} // namespace qmcplusplus

#endif
