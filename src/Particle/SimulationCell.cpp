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


#include "SimulationCell.h"
#include "Lattice/CrystalLattice.h"
#include <type_traits>

namespace qmcplusplus
{

SimulationCell::SimulationCell() = default;

template<typename>
constexpr auto always_false = false;

SimulationCell::SimulationCell(const Lattice& lattice) : lattice_(lattice) { resetLRBox(); }

void SimulationCell::resetLRBox()
{
  if (lattice_.getSuperCellEnum() != SUPERCELL_OPEN)
  {
    lattice_.SetLRCutoffs(lattice_.getRv());
    lrbox_  = lattice_;
    bool changed = false;
    if (lattice_.getSuperCellEnum() == SUPERCELL_SLAB && lattice_.getVacuumScale() != 1.0)
    {
      lrbox_.getR()(2, 0) *= lattice_.getVacuumScale();
      lrbox_.getR()(2, 1) *= lattice_.getVacuumScale();
      lrbox_.getR()(2, 2) *= lattice_.getVacuumScale();
      changed = true;
    }
    else if (lattice_.getSuperCellEnum() == SUPERCELL_WIRE && lattice_.getVacuumScale() != 1.0)
    {
      lrbox_.getR()(1, 0) *= lattice_.getVacuumScale();
      lrbox_.getR()(1, 1) *= lattice_.getVacuumScale();
      lrbox_.getR()(1, 2) *= lattice_.getVacuumScale();
      lrbox_.getR()(2, 0) *= lattice_.getVacuumScale();
      lrbox_.getR()(2, 1) *= lattice_.getVacuumScale();
      lrbox_.getR()(2, 2) *= lattice_.getVacuumScale();
      changed = true;
    }
    lrbox_.reset();
    lrbox_.SetLRCutoffs(lrbox_.getRv());
    lrbox_.printCutoffs(app_log());

    if (changed)
    {
      app_summary() << "  Simulation box changed by vacuum supercell conditions" << std::endl;
      app_log() << "--------------------------------------- " << std::endl;
      lrbox_.print(app_log());
      app_log() << "--------------------------------------- " << std::endl;
    }
    lrbox_.reset();
    lrbox_.SetLRCutoffs(lrbox_.getRv());
    lrbox_.printCutoffs(app_log());

    k_lists_.updateKLists(lrbox_, lrbox_.LR_kc, lrbox_.ndim);
  }
}

} // namespace qmcplusplus
