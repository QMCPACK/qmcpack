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

#include <MinimalParticlePool.h>

namespace qmcplusplus
{
/** This should be the minimal ParticleSetPool for integration tests
 *
 */

ParticleSetPool MinimalParticlePool::make_diamondC_1x1x1(Communicate* c)
{
  ParticleSetPool pp(c);
  parseParticleSetXML(particles_xml, pp);
  return pp;
}

} // namespace qmcplusplus
