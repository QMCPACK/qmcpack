//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_WALKERELEMENTSREF_H
#define QMCPLUSPLUS_WALKERELEMENTSREF_H

#include "Configuration.h"
#include "Particle/Walker.h"

namespace qmcplusplus
{

class TrialWaveFunction;

/** type for returning the walker and its elements from MCPopulation
 *
 *  have no expectations for the validity of the references in this structure past
 *  the context it was returned in. It should not be returned by a call in a
 *  crowd or threaded context.
 * 
 *  @ye-luo's "fat" walker
 *
 *  We need this if we want to "copyFrom" the whole fat walker when it comes off the line
 *  i.e. mpi.  Insuring the "fat" walker is valid at the earliest possible point seems
 *  less likely to end in tears then just calling copyFrom random other places (hopefully)
 *  in time, in order to not access an invalid walker element.
 */
struct WalkerElementsRef
{
  /** to allow use of emplace back
   */
  WalkerElementsRef(Walker<QMCTraits, PtclOnLatticeTraits>& walker_in, ParticleSet& pset_in, TrialWaveFunction& twf_in) : walker(walker_in), pset(pset_in), twf(twf_in) {}
;
  Walker<QMCTraits, PtclOnLatticeTraits>& walker;
  ParticleSet& pset;
  TrialWaveFunction& twf;
};

}

#endif
