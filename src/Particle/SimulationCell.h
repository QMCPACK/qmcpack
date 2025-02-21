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

template<typename REAL>
class SimulationCellT
{
public:
  using FullPrecReal = QMCTraits::FullPrecRealType;
  using Lattice = typename PtclOnLatticeTraits::ParticleLayout;
  using LatticeFullPrec = ParticleLayoutT<FullPrecReal>;
  SimulationCellT();
  // Who uses these, the application uses the default constructor and that fills in lattice_ and full_lattice_
  // from ParticleSetPool
  template<typename REALP, unsigned D>
  SimulationCellT(const CrystalLattice<REALP, D>& lattice);

  const Lattice& getLattice() const { return lattice_; }
  /// Full Precision Lattice, if you doing anything with the lattice wrt kpoints use this!
  const auto& getFullPrecLattice() const { return full_lattice_; }
  const Lattice& getPrimLattice() const { return primative_lattice_; }
  const Lattice& getLRBox() const { return LRBox_; }
  const auto& getFullLRBox() const { return full_lrbox_; }

  void resetLRBox();

  /// access k_lists_ read only
  const KContainer& getKLists() const { return k_lists_; }

  LatticeFullPrec& getFullLattice() { return full_lattice_; }
private:
  ///simulation cell lattice
  Lattice lattice_;
  ///Primative cell lattice
  Lattice primative_lattice_;
  ///long-range box
  Lattice LRBox_;
  ///For kpoint calculations
  LatticeFullPrec full_lattice_;
  LatticeFullPrec full_lrbox_;

  /// K-Vector List.
  KContainer k_lists_;

  friend class ParticleSetPool;
};

using SimulationCell = SimulationCellT<QMCTraits::RealType>;

#ifdef MIXED_PRECISION
extern template class SimulationCellT<float>;
extern template SimulationCellT<float>::SimulationCellT(const CrystalLattice<double, OHMMS_DIM>& lattice);
#else
extern template class SimulationCellT<double>;
extern template SimulationCellT<double>::SimulationCellT(const CrystalLattice<double, OHMMS_DIM>& lattice);
#endif

  
} // namespace qmcplusplus
#endif
