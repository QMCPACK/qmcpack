//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MINIMALPARTICLEPOOL_H
#define QMCPLUSPLUS_MINIMALPARTICLEPOOL_H

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSetPool.h"

namespace qmcplusplus
{
/** This should be the minimal ParticleSetPool for integration tests
 *
 */
class MinimalParticlePool
{
public:
  static void parseParticleSetXML(const char* xml_string, ParticleSetPool& pp);
  static ParticleSetPool make_diamondC_1x1x1(Communicate* comm);
  static ParticleSetPool make_O2_spinor(Communicate* comm);
  static ParticleSetPool make_NiO_a4(Communicate* comm);
  static ParticleSetPool make_H2(Communicate* comm);
};

} // namespace qmcplusplus

#endif
