//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from ParticleSet.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/MoveContext.h"

namespace qmcplusplus
{

MoveContext::MoveContext(int particles)
{
  // This touch might not be necessary
  positions_.resize(particles);
}

void MoveContext::loadWalker(MCPWalker& awalker)
{
  positions_ = awalker.R;
  positions_soa_.copyIn(positions_);
  // in certain cases, full tables must be ready
  // for (int i = 0; i < DistTables.size(); i++)
  //     if (DistTables[i]->DTType == DT_AOS || DistTables[i]->Need_full_table_loadWalker)
  //       DistTables[i]->evaluate(*this);
  //   //computed so that other objects can use them, e.g., kSpaceJastrow
  //   if (SK && SK->DoUpdate)
  //     SK->UpdateAllPart(*this);

  // activePtcl = -1;

}

}
