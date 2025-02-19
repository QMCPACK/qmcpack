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

ParticleLayoutT<QMCTraits::FullPrecRealType> makeFullPrecParticleLayout(xmlNodePtr cur);


class LatticeParser
{
  ParticleLayoutT<QMCTraits::FullPrecRealType>& full_prec_ref_;
  using ParticleLayout = ParticleLayoutT<QMCTraits::RealType>;
  ParticleLayout& ref_;
public:
  template<typename T1, unsigned D>
  LatticeParser(CrystalLattice<T1,D>& full_lat) : full_prec_ref_(full_lat), ref_(full_lat) {
    static_assert(std::is_same_v<T1, QMCTraits::FullPrecRealType>, "This should never get called with a reduced precision lattice");
  }
  template<typename T1, typename T2, unsigned D>
  LatticeParser(CrystalLattice<T1,D>& full_lat, CrystalLattice<T2,D>& lat ) : full_prec_ref_(full_lat), ref_(lat) {}
  bool put(xmlNodePtr cur);
};

class LatticeXMLWriter
{
  template<typename T>
  using ParticleLayoutT = typename PtclOnLatticeTraitsT<T>::ParticleLayout;
  using ParticleLayout = ParticleLayoutT<QMCTraits::FullPrecRealType>;
  const ParticleLayout& ref_;

public:
  LatticeXMLWriter(const ParticleLayout& lat) : ref_(lat) {}
  bool get(std::ostream&) const;
  xmlNodePtr createNode();
};

  
} // namespace qmcplusplus
#endif
