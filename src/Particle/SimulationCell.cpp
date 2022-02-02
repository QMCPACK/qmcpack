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


#include "SimulationCell.h"

namespace qmcplusplus
{

SimulationCell::SimulationCell() = default;

SimulationCell::SimulationCell(const Lattice& lattice)
    : lattice_(lattice)
{
  resetLRBox();
}

void SimulationCell::resetLRBox()
{
  if (lattice_.SuperCellEnum != SUPERCELL_OPEN)
  {
    lattice_.SetLRCutoffs(lattice_.Rv);
    LRBox_        = lattice_;
    bool changed = false;
    if (lattice_.SuperCellEnum == SUPERCELL_SLAB && lattice_.VacuumScale != 1.0)
    {
      LRBox_.R(2, 0) *= lattice_.VacuumScale;
      LRBox_.R(2, 1) *= lattice_.VacuumScale;
      LRBox_.R(2, 2) *= lattice_.VacuumScale;
      changed = true;
    }
    else if (lattice_.SuperCellEnum == SUPERCELL_WIRE && lattice_.VacuumScale != 1.0)
    {
      LRBox_.R(1, 0) *= lattice_.VacuumScale;
      LRBox_.R(1, 1) *= lattice_.VacuumScale;
      LRBox_.R(1, 2) *= lattice_.VacuumScale;
      LRBox_.R(2, 0) *= lattice_.VacuumScale;
      LRBox_.R(2, 1) *= lattice_.VacuumScale;
      LRBox_.R(2, 2) *= lattice_.VacuumScale;
      changed = true;
    }
    LRBox_.reset();
    LRBox_.SetLRCutoffs(LRBox_.Rv);
    LRBox_.printCutoffs(app_log());

    if (changed)
    {
      app_summary() << "  Simulation box changed by vacuum supercell conditions" << std::endl;
      app_log() << "--------------------------------------- " << std::endl;
      LRBox_.print(app_log());
      app_log() << "--------------------------------------- " << std::endl;
    }

    k_lists_.updateKLists(LRBox_, LRBox_.LR_kc);
  }
}
}
