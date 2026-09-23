//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#include "L2Diffusion.h"

#include <algorithm>
#include <cassert>

#include "Containers/OhmmsPETE/TensorOps.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{
void L2Diffusion::Workspace::resize(size_t num_walkers)
{
  if (move_valid.size() == num_walkers)
  {
    assert(zero_displacements.positions.size() == num_walkers);
    assert(reject_all_intermediate.size() == num_walkers);
    assert(diffusion_tensors.size() == num_walkers);
    assert(drift_corrections.size() == num_walkers);
    return;
  }

  // MCCoords positions and vector<bool> value-initialize new entries to zero
  // and false, respectively. Both buffers are passed as const input afterward.
  zero_displacements.positions.resize(num_walkers);
  reject_all_intermediate.resize(num_walkers, false);
  move_valid.resize(num_walkers);
  diffusion_tensors.resize(num_walkers);
  drift_corrections.resize(num_walkers);
}

L2Diffusion::PosType L2Diffusion::computeScaledDrift(RealType tauovermass,
                                                     const GradType& gradient,
                                                     const TensorType& diffusion_tensor,
                                                     PosType drift_correction)
{
  PosType drift;
  getScaledDriftL2(tauovermass, gradient, diffusion_tensor, drift_correction, drift);
  return drift;
}

void L2Diffusion::addDiffusion(const TensorType& diffusion_tensor,
                               PosType& gaussian_displacement,
                               PosType& proposed_displacement)
{
  gaussian_displacement = dot(cholesky(diffusion_tensor), gaussian_displacement);
  proposed_displacement += gaussian_displacement;
}

void L2Diffusion::prepareMove(const TauParams<RealType, CoordsType::POS>& taus,
                              const TWFGrads<CoordsType::POS>& grads_now,
                              int iat,
                              const PSdispatcher& ps_dispatcher,
                              const Hdispatcher& ham_dispatcher,
                              const RefVectorWithLeader<ParticleSet>& walker_elecs,
                              const RefVectorWithLeader<QMCHamiltonian>& walker_hamiltonians,
                              MCCoords<CoordsType::POS>& gaussian_displacements,
                              MCCoords<CoordsType::POS>& proposed_displacements,
                              std::vector<RealType>& log_gf,
                              std::vector<bool>& are_valid,
                              Workspace& workspace) const
{
  const size_t num_walkers = walker_elecs.size();
  assert(walker_hamiltonians.size() == num_walkers);
  assert(gaussian_displacements.positions.size() == num_walkers);
  assert(proposed_displacements.positions.size() == num_walkers);
  assert(grads_now.grads_positions.size() == num_walkers);
  assert(log_gf.size() == num_walkers);
  assert(are_valid.size() == num_walkers);

  workspace.resize(num_walkers);

  for (auto& displacement : gaussian_displacements.positions)
    displacement *= taus.sqrttau;
  std::transform(gaussian_displacements.positions.begin(), gaussian_displacements.positions.end(), log_gf.begin(),
                 [halfovertau = taus.oneover2tau](const PosType& displacement) {
                   return -halfovertau * dot(displacement, displacement);
                 });

  // L2Potential evaluates D and K from temporary distance-table data. A zero
  // move stages the current position in those buffers.
  ps_dispatcher.flex_makeMove(walker_elecs, iat, workspace.zero_displacements, are_valid);
  for (size_t iw = 0; iw < num_walkers; ++iw)
    workspace.move_valid[iw] = are_valid[iw];
  ham_dispatcher.flex_computeL2DK(walker_hamiltonians, walker_elecs, iat, workspace.diffusion_tensors,
                                  workspace.drift_corrections);
  ps_dispatcher.flex_accept_rejectMove<CoordsType::POS>(walker_elecs, iat, workspace.reject_all_intermediate);

  for (size_t iw = 0; iw < num_walkers; ++iw)
    proposed_displacements.positions[iw] =
        computeScaledDrift(taus.tauovermass, grads_now.grads_positions[iw], workspace.diffusion_tensors[iw],
                           workspace.drift_corrections[iw]);

  // Evaluate the diffusion tensor at the drifted position.
  ps_dispatcher.flex_makeMove(walker_elecs, iat, proposed_displacements, are_valid);
  for (size_t iw = 0; iw < num_walkers; ++iw)
    workspace.move_valid[iw] = workspace.move_valid[iw] && are_valid[iw];
  ham_dispatcher.flex_computeL2D(walker_hamiltonians, walker_elecs, iat, workspace.diffusion_tensors);
  ps_dispatcher.flex_accept_rejectMove<CoordsType::POS>(walker_elecs, iat, workspace.reject_all_intermediate);

  for (size_t iw = 0; iw < num_walkers; ++iw)
    addDiffusion(workspace.diffusion_tensors[iw], gaussian_displacements.positions[iw],
                 proposed_displacements.positions[iw]);
}

void L2Diffusion::applyMoveValidity(std::vector<bool>& are_valid, const Workspace& workspace) const
{
  assert(are_valid.size() == workspace.move_valid.size());
  for (size_t iw = 0; iw < are_valid.size(); ++iw)
    are_valid[iw] = are_valid[iw] && workspace.move_valid[iw];
}
} // namespace qmcplusplus
