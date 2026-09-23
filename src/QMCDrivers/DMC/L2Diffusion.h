//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_L2DIFFUSION_H
#define QMCPLUSPLUS_L2DIFFUSION_H

#include "Configuration.h"
#include "Particle/MCCoords.hpp"
#include "Particle/PSdispatcher.h"
#include "QMCHamiltonians/Hdispatcher.h"
#include "QMCDrivers/TauParams.hpp"
#include "QMCWaveFunctions/TWFGrads.hpp"

namespace qmcplusplus
{
/** Encapsulates construction of particle proposals for an L2 pseudopotential.
 *
 * The object is immutable and may be shared by crowd tasks. Mutable scratch data
 * is kept in one Workspace per crowd.
 */
class L2Diffusion : public QMCTraits
{
public:
  using TensorType = QMCHamiltonian::TensorType;

  struct Workspace
  {
    explicit Workspace(size_t num_walkers = 0) { resize(num_walkers); }

    void resize(size_t num_walkers);

    /// Immutable false-valued decisions used to reject temporary moves.
    std::vector<bool> reject_all_intermediate;
    /// Output scratch; each entry is overwritten before it is read.
    std::vector<bool> move_valid;
    /// Output scratch; the Hamiltonian overwrites each entry before it is read.
    std::vector<TensorType> diffusion_tensors;
    /// Output scratch; the Hamiltonian overwrites each entry before it is read.
    std::vector<PosType> drift_corrections;
  };

  /** Build the L2 drift-diffusion displacement for one particle across a crowd. */
  static void prepareMove(const TauParams<RealType, CoordsType::POS>& taus,
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
                   Workspace& workspace);

  /** Combine final-move validity with validity of the temporary L2 moves. */
  static void applyMoveValidity(std::vector<bool>& are_valid, const Workspace& workspace);

  /** Compute the scaled drift used by the L2 propagator. */
  static PosType computeScaledDrift(RealType tauovermass,
                                    const GradType& gradient,
                                    const TensorType& diffusion_tensor,
                                    PosType drift_correction);

  /** Transform a Gaussian displacement by the L2 diffusion tensor and add it to the drift. */
  static void addDiffusion(const TensorType& diffusion_tensor,
                           PosType& gaussian_displacement,
                           PosType& proposed_displacement);
};
} // namespace qmcplusplus

#endif
