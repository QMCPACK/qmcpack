//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SIMULATIONCELL_H
#define QMCPLUSPLUS_SIMULATIONCELL_H

#include "Configuration.h"
#include "LongRange/KContainer.h"

namespace qmcplusplus
{
class ParticleSetPool;

class SimulationCell
{
public:
  using Lattice = PtclOnLatticeTraits::ParticleLayout;

  SimulationCell();
  SimulationCell(const Lattice& lattice);

  const Lattice& getLattice() const { return lattice_; }
  const Lattice& getPrimLattice() const { return primative_lattice_; }
  const Lattice& getLRBox() const { return LRBox_; }

  void resetLRBox();

  /// access k_lists_ read only
  const KContainer& getKLists() const { return k_lists_; }

private:
  ///simulation cell lattice
  Lattice lattice_;
  ///Primative cell lattice
  Lattice primative_lattice_;
  ///long-range box
  Lattice LRBox_;

  /// K-Vector List.
  KContainer k_lists_;

  friend class ParticleSetPool;
};
} // namespace qmcplusplus
#endif
