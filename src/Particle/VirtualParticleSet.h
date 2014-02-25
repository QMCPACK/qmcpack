//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by QMCPACK developers
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

    vector<ValueType> ratios;
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
/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5985 $   $Date: 2013-09-27 17:08:27 -0400 (Fri, 27 Sep 2013) $
 * $Id: ParticleSet.h 5985 2013-09-27 21:08:27Z jtkrogel $
 ***************************************************************************/
