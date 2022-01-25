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

#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

using Lattice = ParticleSet::ParticleLayout;

/** create super lattice */
Lattice createSuperLattice(const Lattice& in, const Tensor<int, OHMMS_DIM>& tmat);
/** expand a particleset */
void expandSuperCell(const ParticleSet& in, const Tensor<int, OHMMS_DIM>& tmat, ParticleSet& out);

} // namespace qmcplusplus
#endif
