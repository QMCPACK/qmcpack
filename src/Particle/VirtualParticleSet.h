//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
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
  class VirtualParticleSet:  public ParticleSet
  {
    /// ParticleSet this object refers to
    const ParticleSet* myPtcl;
    /** initialize minimum data
     *
     * Create DistTables of AB type
     */
    void init_minimum(int np);
    public:

    std::vector<ValueType> ratios;
    /** constructor 
     * @param p ParticleSet whose virtual moves are handled by this object
     * @param nptcl number of virtual particles
     */
    VirtualParticleSet(ParticleSet* p, int nptcl=0);

    ~VirtualParticleSet();

    /** move the iat-th particle of myPtcl by multiple displacements
     *
     * DistTables[0]=dist(myPtcl,this)
     * DistTables[other]=dist(other,this)
     */
    void makeMoves(int iat, const ParticlePos_t& displ);

    void validate(int iel, int k);

    inline const DistanceTableData* getVirtualTable(int i) const
    {
      return DistTables[i];
    }

    inline const DistanceTableData* getRealTable(int i) const
    {
      return myPtcl->DistTables[i];
    }

    void reset(const ParticleSet* p);
  };
}
#endif
