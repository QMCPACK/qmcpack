//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <vector>
#include <complex>
#include "Configuration.h"
#include "ParticleSet.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus
{
/** evalaute bare Coulomb in 3D using LRHandlerTemp
 */
TEST_CASE("StructFact", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Approx(Lattice.Volume) == 125);
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Approx(Lattice.LR_rc) == 2.5);
  REQUIRE(Approx(Lattice.LR_kc) == 12);

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.SK.getKLists()

  SpeciesSet& tspecies = ref.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  int downIdx          = tspecies.addSpecies("d");

  ref.create({3,1});
  ref.R[0] = {0.0, 1.0, 2.0};
  ref.R[1] = {1.0, 0.2, 3.0};
  ref.R[2] = {0.3, 4.0, 1.4};
  ref.R[3] = {3.2, 4.7, 0.7};

  StructFact sk(tspecies.size(), ref.getTotalNum(), ref.getLRBox(), 50);
  REQUIRE(sk.getKLists().numk == 263786);
  sk.updateAllPart(ref);

  std::vector<std::complex<double>> rhok_sum_ref{-125.80618630936, 68.199075127271};

  for (int i = 0; i < ref.groups(); i++)
  {
    std::complex<QMCTraits::RealType> rhok_sum, rhok_even_sum;
    for (int ik = 0; ik < sk.getKLists().numk; ik++)
      rhok_sum += std::complex<QMCTraits::RealType>(sk.rhok_r[i][ik], sk.rhok_i[i][ik]);

    //std::cout << std::setprecision(14) << rhok_sum << std::endl;
    CHECK(ComplexApprox(rhok_sum).epsilon(5e-5) == rhok_sum_ref[i]);
  }
}

} // namespace qmcplusplus
