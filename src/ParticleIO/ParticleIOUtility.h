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


#ifndef OHMMS_PARTICLE_INPUTOUTPUT_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
/** expand a particleset including lattice */
void expandSuperCell(ParticleSet& ref, Tensor<int, OHMMS_DIM>& tmat);

} // namespace qmcplusplus
#endif
