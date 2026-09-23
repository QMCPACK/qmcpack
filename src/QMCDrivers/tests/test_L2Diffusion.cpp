//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "QMCDrivers/DMC/L2Diffusion.h"

namespace qmcplusplus
{
TEST_CASE("L2Diffusion workspace and proposal math", "[drivers][l2]")
{
  using RealType = QMCTraits::RealType;
  using PosType  = QMCTraits::PosType;

  L2Diffusion diffusion;
  L2Diffusion::Workspace workspace(2);
  REQUIRE(workspace.zero_displacements.positions.size() == 2);
  REQUIRE(workspace.reject_all_intermediate.size() == 2);
  REQUIRE(workspace.move_valid.size() == 2);
  REQUIRE(workspace.diffusion_tensors.size() == 2);
  REQUIRE(workspace.drift_corrections.size() == 2);
  for (const PosType& displacement : workspace.zero_displacements.positions)
    CHECK(displacement == PosType(0.0));
  CHECK(workspace.reject_all_intermediate == std::vector<bool>{false, false});

  workspace.resize(3);
  REQUIRE(workspace.zero_displacements.positions.size() == 3);
  REQUIRE(workspace.reject_all_intermediate.size() == 3);
  REQUIRE(workspace.move_valid.size() == 3);
  REQUIRE(workspace.diffusion_tensors.size() == 3);
  REQUIRE(workspace.drift_corrections.size() == 3);
  for (const PosType& displacement : workspace.zero_displacements.positions)
    CHECK(displacement == PosType(0.0));
  CHECK(workspace.reject_all_intermediate == std::vector<bool>{false, false, false});

  workspace.resize(2);

  L2Diffusion::TensorType drift_tensor;
  drift_tensor = 0.0;
  drift_tensor(0, 0) = 2.0;
  drift_tensor(1, 1) = 3.0;
  drift_tensor(2, 2) = 4.0;
  const QMCTraits::GradType gradient{1.0, -2.0, 0.5};
  const PosType drift_correction{0.5, 0.25, -1.0};
  constexpr RealType tau = 0.2;

  const PosType bare_drift{1.5, -6.25, 3.0};
  const RealType bare_norm = dot(bare_drift, bare_drift);
  const RealType scale     = (-1.0 + std::sqrt(1.0 + 2.0 * tau * bare_norm)) / bare_norm;
  const PosType drift      = L2Diffusion::computeScaledDrift(tau, gradient, drift_tensor, drift_correction);
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
    CHECK(drift[idim] == Approx(bare_drift[idim] * scale));

  L2Diffusion::TensorType diffusion_tensor;
  diffusion_tensor = 0.0;
  diffusion_tensor(0, 0) = 4.0;
  diffusion_tensor(1, 1) = 9.0;
  diffusion_tensor(2, 2) = 16.0;
  PosType gaussian{1.0, 2.0, -1.0};
  PosType proposal{0.25, 0.5, 0.75};
  L2Diffusion::addDiffusion(diffusion_tensor, gaussian, proposal);
  CHECK(gaussian[0] == Approx(2.0));
  CHECK(gaussian[1] == Approx(6.0));
  CHECK(gaussian[2] == Approx(-4.0));
  CHECK(proposal[0] == Approx(2.25));
  CHECK(proposal[1] == Approx(6.5));
  CHECK(proposal[2] == Approx(-3.25));

  workspace.move_valid = {true, false};
  std::vector<bool> final_move_valid{false, true};
  diffusion.applyMoveValidity(final_move_valid, workspace);
  CHECK(final_move_valid == std::vector<bool>{false, false});
}
} // namespace qmcplusplus
