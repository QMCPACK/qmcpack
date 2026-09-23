//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "QMCHamiltonians/L2Potential.h"

namespace qmcplusplus
{
TEST_CASE("L2Potential diffusion and drift correction", "[hamiltonian][l2]")
{
  using RealType   = QMCTraits::RealType;
  using PosType    = QMCTraits::PosType;
  using TensorType = QMCTraits::TensorType;

  Lattice lattice;
  lattice.R.diagonal(10.0);
  lattice.reset();
  const SimulationCell simulation_cell(lattice);

  ParticleSet ions(simulation_cell);
  ions.setName("ion");
  ions.create({1});
  ions.getSpeciesSet().addSpecies("X");
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.update();

  ParticleSet electrons(simulation_cell);
  electrons.setName("e");
  electrons.create({1});
  electrons.R[0] = {1.0, 0.0, 0.0};

  L2Potential l2_potential(ions, electrons);
  electrons.update();

  auto grid = std::make_unique<LinearGrid<RealType>>();
  grid->set(0.0, 2.0, 5);
  std::vector<RealType> values(5, 0.5);
  auto radial_potential = std::make_unique<L2RadialPotential>();
  radial_potential->vL2 =
      std::make_unique<L2RadialPotential::RadialPotentialType>(std::move(grid), values);
  radial_potential->vL2->spline();
  radial_potential->rcut = 2.0;
  l2_potential.add(0, std::move(radial_potential));

  TensorType diffusion_tensor;
  PosType drift_correction;
  l2_potential.evaluateDK(electrons, 0, diffusion_tensor, drift_correction);

  CHECK(diffusion_tensor(0, 0) == Approx(1.0));
  CHECK(diffusion_tensor(1, 1) == Approx(2.0));
  CHECK(diffusion_tensor(2, 2) == Approx(2.0));
  CHECK(diffusion_tensor(0, 1) == Approx(0.0));
  CHECK(diffusion_tensor(0, 2) == Approx(0.0));
  CHECK(diffusion_tensor(1, 2) == Approx(0.0));
  CHECK(drift_correction[0] == Approx(1.0));
  CHECK(drift_correction[1] == Approx(0.0));
  CHECK(drift_correction[2] == Approx(0.0));

  TensorType diffusion_only;
  electrons.makeMove(0, PosType{0.0, 0.0, 0.0});
  l2_potential.evaluateD(electrons, 0, diffusion_only);
  for (int i = 0; i < OHMMS_DIM; ++i)
    for (int j = 0; j < OHMMS_DIM; ++j)
      CHECK(diffusion_only(i, j) == Approx(diffusion_tensor(i, j)));
  electrons.rejectMove(0);

  // Outside the cutoff, L2 contributes neither diffusion nor drift correction.
  electrons.makeMove(0, PosType{2.0, 0.0, 0.0});
  electrons.acceptMove(0);
  l2_potential.evaluateDK(electrons, 0, diffusion_tensor, drift_correction);
  for (int i = 0; i < OHMMS_DIM; ++i)
    for (int j = 0; j < OHMMS_DIM; ++j)
      CHECK(diffusion_tensor(i, j) == Approx(i == j ? 1.0 : 0.0));
  CHECK(drift_correction == PosType(0.0));
}
} // namespace qmcplusplus
