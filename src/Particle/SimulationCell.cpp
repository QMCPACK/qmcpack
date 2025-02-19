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


#include "SimulationCell.h"
#include "Lattice/CrystalLattice.h"
#include <type_traits>

namespace qmcplusplus
{

template<typename REAL>
SimulationCellT<REAL>::SimulationCellT() = default;

  template<typename>
  constexpr auto always_false = false;
  
template<typename REAL>
template<typename REALP, unsigned D>
SimulationCellT<REAL>::SimulationCellT(const CrystalLattice<REALP, D>& lattice)
{
  static_assert(std::is_floating_point_v<REAL> && std::is_floating_point_v<REALP> );
  if constexpr (!std::is_same_v<REALP, QMCTraits::FullPrecRealType>)
    static_assert(always_false<REALP>, "You should never be trying to use a reduced precision lattice to construct a simulation cell");
  full_lattice_ = lattice;
  resetLRBox();
}
  
template<typename REAL>
void SimulationCellT<REAL>::resetLRBox()
{
  // We don't want to use the reduced precision lattice for any of this
  if (full_lattice_.getSuperCellEnum() != SUPERCELL_OPEN)
  {
    full_lattice_.SetLRCutoffs(full_lattice_.getRv());
    full_lrbox_       = full_lattice_;
    bool changed = false;
    if (full_lattice_.getSuperCellEnum() == SUPERCELL_SLAB && full_lattice_.getVacuumScale() != 1.0)
    {
      full_lrbox_.getR()(2, 0) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(2, 1) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(2, 2) *= full_lattice_.getVacuumScale();
      changed = true;
    }
    else if (full_lattice_.getSuperCellEnum() == SUPERCELL_WIRE && full_lattice_.getVacuumScale() != 1.0)
    {
      full_lrbox_.getR()(1, 0) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(1, 1) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(1, 2) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(2, 0) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(2, 1) *= full_lattice_.getVacuumScale();
      full_lrbox_.getR()(2, 2) *= full_lattice_.getVacuumScale();
      changed = true;
    }
    full_lrbox_.reset();
    full_lrbox_.SetLRCutoffs(full_lrbox_.getRv());
    full_lrbox_.printCutoffs(app_log());

    if (changed)
    {
      app_summary() << "  Simulation box changed by vacuum supercell conditions" << std::endl;
      app_log() << "--------------------------------------- " << std::endl;
      full_lrbox_.print(app_log());
      app_log() << "--------------------------------------- " << std::endl;
    }
    full_lrbox_.reset();
    full_lrbox_.SetLRCutoffs(full_lrbox_.getRv());
    full_lrbox_.printCutoffs(app_log());
    
    k_lists_.updateKLists<FullPrecReal>(full_lrbox_, full_lrbox_.LR_kc, full_lrbox_.ndim);
  }
  LRBox_ = full_lrbox_;
  lattice_ = full_lattice_;
}

#ifdef MIXED_PRECISION
template class SimulationCellT<float>;
template SimulationCellT<float>::SimulationCellT(const CrystalLattice<double, OHMMS_DIM>& lattice);
#else
template class SimulationCellT<double>;
template SimulationCellT<double>::SimulationCellT(const CrystalLattice<double, OHMMS_DIM>& lattice);
#endif
} // namespace qmcplusplus
