//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NEIGHBORLIST_H
#define QMCPLUSPLUS_NEIGHBORLIST_H

#include <vector>
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
class NeighborLists
{
protected:
  /// neighbor particle IDs of all the reference particles
  std::vector<std::vector<int>> NeighborIDs;

public:
  /// constructor
  NeighborLists(const ParticleSet& source) : NeighborIDs(source.getTotalNum()) {}

  /// get the neighbor list of the source particle
  std::vector<int>& getNeighborList(int source) { return NeighborIDs[source]; }
  const std::vector<int>& getNeighborList(int source) const { return NeighborIDs[source]; }
};

} // namespace qmcplusplus
#endif
