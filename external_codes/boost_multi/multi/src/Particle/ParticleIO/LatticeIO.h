//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_PARTICLELAYOUT_INPUTOUTPUT_UTILITY_H
#define OHMMS_PARTICLELAYOUT_INPUTOUTPUT_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Configuration.h"

namespace qmcplusplus
{
class LatticeParser
{
  using ParticleLayout = PtclOnLatticeTraits::ParticleLayout;
  ParticleLayout& ref_;

public:
  LatticeParser(ParticleLayout& lat) : ref_(lat) {}
  bool put(xmlNodePtr cur);
};


class LatticeXMLWriter
{
  using ParticleLayout = PtclOnLatticeTraits::ParticleLayout;
  const ParticleLayout& ref_;

public:
  LatticeXMLWriter(const ParticleLayout& lat) : ref_(lat) {}
  bool get(std::ostream&) const;
  xmlNodePtr createNode();
};


} // namespace qmcplusplus
#endif
