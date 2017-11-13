//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file VirtualParticleSet.h
 * A proxy class to the quantum ParticleSet
 */
#ifndef QMCPLUSPLUS_VIRTUAL_PARTICLESET_H
#define QMCPLUSPLUS_VIRTUAL_PARTICLESET_H

#include <Configuration.h>
#include <Particle/ParticleSet.h>

namespace qmcplusplus
{

  /** Introduced to handle virtual moves and ratio computations, e.g. for non-local PP evaluations.
   */
  class VirtualParticleSet: public ParticleSet
  {
    public:
    ///reference particle ID
    Index_t refID;

    /// ParticleSet this object refers to
    const ParticleSet& refPtcl;

    /** constructor 
     * @param p ParticleSet whose virtual moves are handled by this object
     * @param nptcl number of virtual particles
     */
    VirtualParticleSet(const ParticleSet& p, int nptcl);

    /// move virtual particles to new postions and update distance tables
    void makeMoves(int iat, const ParticlePos_t& vitualPos);
  };
}
#endif
