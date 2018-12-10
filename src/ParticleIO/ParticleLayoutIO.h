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

  typedef PtclOnLatticeTraits::ParticleLayout_t ParticleLayout_t;
  ParticleLayout_t& ref_;

public:

  LatticeParser(ParticleLayout_t& lat): ref_(lat) { }
  bool put(xmlNodePtr cur);
};


class LatticeXMLWriter
{

  typedef PtclOnLatticeTraits::ParticleLayout_t ParticleLayout_t;
  ParticleLayout_t& ref_;
public:

  LatticeXMLWriter(ParticleLayout_t& lat): ref_(lat) { }
  bool get(std::ostream& ) const;
  xmlNodePtr createNode();
};


}
#endif

