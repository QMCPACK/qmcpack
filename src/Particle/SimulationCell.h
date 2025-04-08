//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
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
  using FullPrecReal = QMCTraits::FullPrecRealType;
  SimulationCell();
  SimulationCell(const Lattice& lattice);

  const Lattice& getLattice() const { return lattice_; }
  const Lattice& getPrimLattice() const { return primitive_lattice_; }
  const Lattice& getLRBox() const { return lrbox_; }
  const KContainer& getKLists() const { return k_lists_; }

  Lattice& getLattice() { return lattice_; }

  void resetLRBox();
private:
  ///simulation cell lattice
  Lattice lattice_;
  ///Primative cell lattice
  Lattice primitive_lattice_;
  ///long-range box
  Lattice lrbox_;

  /// K-Vector List.
  KContainer k_lists_;

  friend class ParticleSetPool;
};

} // namespace qmcplusplus
#endif
