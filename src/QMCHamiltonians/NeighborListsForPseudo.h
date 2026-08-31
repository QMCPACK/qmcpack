//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NEIGHBORLIST_H
#define QMCPLUSPLUS_NEIGHBORLIST_H

#include <vector>
#include "NonLocalECPComponent.h"
#include <DistanceTable.h>

namespace qmcplusplus
{
class NeighborLists
{
  /// neighboring particle IDs of all the reference particles
  std::vector<std::vector<int>> neighborIDs_;

public:
  /** constructor
   * @param num_ref number of reference particles
   */
  NeighborLists(size_t num_ref) : neighborIDs_(num_ref) {}

  /// get the number of neignbor lists / reference particles.
  size_t size() const { return neighborIDs_.size(); }
  /// get the neighbor list of the source particle
  std::vector<int>& getNeighborList(int source) { return neighborIDs_[source]; }
  const std::vector<int>& getNeighborList(int source) const { return neighborIDs_[source]; }
};

class NeighborListsForPseudo
{
  ///neighborlist of electrons
  NeighborLists elec_neighbor_ions_;
  ///neighborlist of ions
  NeighborLists ion_neighbor_elecs_;
  ///the set of local-potentials (one for each ion)
  const std::vector<NonLocalECPComponent*>& PP;

public:
  NeighborListsForPseudo(size_t num_elecs, size_t num_ions, const std::vector<NonLocalECPComponent*>& pp);

  /// get the neighboring ion list of a given electron
  const std::vector<int>& getNeighboringIons(int jel) const { return elec_neighbor_ions_.getNeighborList(jel); }

  /** mark all the electrons affected by T-moves and update elec_neighbor_ions_ and ion_neighbor_elecs_
   * @param myTable electron ion distance table
   * @param iel reference electron
   * Note this function should be called before acceptMove for a Tmove
   */
  void markAffectedElecs(const DistanceTableAB& myTable, int iel, std::vector<bool>& elecTMAffected);

  /// clear all the electron and ion neighbor lists
  void clear();

  /** add electron ion pair to the neighbor lists
   * @param jel electron index
   * @param iat ion index
   */
  void addElecIonPair(int jel, int iat);
};

} // namespace qmcplusplus
#endif
