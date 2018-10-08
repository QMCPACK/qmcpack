//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   refactored from WaveFunctionComponent.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WAVEFUNCTIONCOMPONENTTYPEALIASES_H
#define QMCPLUSPLUS_WAVEFUNCTIONCOMPONENTTYPEALIASES_H

#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"

namespace qmcplusplus
{

/*! Collects common types for WaveFunctionComponent consumers and children
 *  When MI and the idea of gaining typdefs and typealiases through
 *  inheritance ambiguous name resolution is soon to follow.
 *  This prevents that and can help make more explicit where
 *  types actually originate.
 */

struct WaveFunctionComponentTypeAliases
{
  typedef ParticleSet::Walker_t     Walker_t;
  typedef Walker_t::WFBuffer_t      WFBufferType;
  typedef Walker_t::Buffer_t        BufferType;
};

}
#endif
