//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <vector>
#include <iostream>
#include "ParticleSet.h"
#include "Lattice/ParticleBConds3DSoa.h"
#include "SoaDistanceTableAA.h"

namespace qmcplusplus
{
TEST_CASE("SoaDistanceTableAA compute_size", "[distance_table]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);

  elec.setName("e");
  elec.create({6, 4});
  
  // using open BC
  SoaDistanceTableAA<OHMMS_PRECISION, OHMMS_DIM, SUPERCELL_OPEN + SOA_OFFSET> dt_ee(elec);

  const size_t Alignment = getAlignment<OHMMS_PRECISION>();

  // create reference values
  assert(elec.getTotalNum() == 10);
  std::vector<int> ref_results(10);
  if (Alignment == 4)
    ref_results = {0, 4, 8, 12, 16, 24, 32, 40, 48, 60};
  else if (Alignment == 8)
    ref_results = {0, 8, 16, 24, 32, 40, 48, 56, 64, 80};

  // run checks
  if (Alignment == 4 || Alignment == 8)
  {
    std::cout << "testing Alignment = " << Alignment << std::endl;
    for (int i = 0; i < ref_results.size(); i++)
      CHECK(dt_ee.compute_size(i) == ref_results[i]);
  }
}
} // namespace qmcplusplus
