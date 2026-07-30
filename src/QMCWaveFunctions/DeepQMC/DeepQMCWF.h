//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/** @file
 * @brief Declaration of a wavefunction component backed by DeepQMC neural-network inference.
 */
#ifndef QMCPLUSPLUS_DEEPQMCWF_H
#define QMCPLUSPLUS_DEEPQMCWF_H

#include <memory>
#include <string>
#include <vector>

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCBridge.h"

namespace qmcplusplus
{

/** Wavefunction component that evaluates a DeepQMC neural-network ansatz.
 *
 * DeepQMCWF delegates log-value, gradient, and laplacian evaluation to a
 * DeepQMCBridge implementation.  The component stores a reference to the ion
 * ParticleSet, gathers electron coordinates from each target ParticleSet, and
 * batches multi-walker calls through the shared bridge object.
 *
 * Particle-by-particle ratio and gradient methods evaluate the proposed
 * electron position through the bridge and cache the proposed log value until
 * acceptMove() or restore() completes the move.
 */
class DeepQMCWF : public WaveFunctionComponent
{
  /// Private tag selecting the shared-bridge constructor used only by makeClone().
  struct SharedBridgeTag
  {};

public:
  /** Constructor transferring ownership of the DeepQMC bridge.
   * @param name component name from the input wavefunction node
   * @param ions source ion ParticleSet used to provide nuclear coordinates
   * @param bridge bridge object used for DeepQMC inference
   * @param mol_idx molecule index passed through to the DeepQMC model
   */
  DeepQMCWF(std::string name, const ParticleSet& ions, std::unique_ptr<const DeepQMCBridge> bridge, int mol_idx);

  /** Constructor sharing an existing DeepQMC bridge for makeClone().
   *
   * The leading SharedBridgeTag makes this overload unambiguous when callers pass
   * a std::unique_ptr to a concrete DeepQMCBridge implementation and documents
   * that shared-bridge construction is not part of the normal external API.  The
   * tag type is private, so outside code cannot name or create it, but the
   * constructor itself is public so std::make_unique can instantiate it from
   * makeClone().
   *
   * @param name component name from the input wavefunction node
   * @param ions source ion ParticleSet used to provide nuclear coordinates
   * @param bridge shared bridge object used for DeepQMC inference
   * @param mol_idx molecule index passed through to the DeepQMC model
   */
  DeepQMCWF(SharedBridgeTag,
            std::string name,
            const ParticleSet& ions,
            std::shared_ptr<const DeepQMCBridge> bridge,
            int mol_idx);

  /// Return the concrete wavefunction component class name.
  std::string getClassName() const override { return "DeepQMCWF"; }

  LogValue evaluateLog(const ParticleSet& P,
                       ParticleSet::ParticleGradient& G,
                       ParticleSet::ParticleLaplacian& L) override;

  void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                      const RefVectorWithLeader<ParticleSet>& p_list,
                      const RefVector<ParticleSet::ParticleGradient>& G_list,
                      const RefVector<ParticleSet::ParticleLaplacian>& L_list) const override;

  void mw_prepareGroup(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                       const RefVectorWithLeader<ParticleSet>& p_list,
                       int ig) const override;
  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   int iat,
                   std::vector<GradType>& grad_now) const override;
  void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios) const override;
  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValue>& ratios,
                    std::vector<GradType>& grad_new) const override;
  void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            const std::vector<bool>& isAccepted,
                            bool safe_to_delay) const override;
  void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const override;
  void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                     const RefVectorWithLeader<ParticleSet>& p_list,
                     const RefVector<ParticleSet::ParticleGradient>& G_list,
                     const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                     bool fromscratch) const override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;
  void restore(int iat) override;
  PsiValue ratio(ParticleSet& P, int iat) override;
  GradType evalGrad(ParticleSet& P, int iat) override;
  PsiValue ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override {}
  LogValue updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override {}

  void evaluateDerivatives(ParticleSet& P,
                           const OptVariables& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tpq) const override;

  /// Return the ion ParticleSet supplying nuclear coordinates.
  const ParticleSet& getIons() const { return ions_; }
  /// Return the bridge used for DeepQMC inference.
  const DeepQMCBridge& getBridge() const { return *bridge_; }
  /// Return the molecule index passed to DeepQMC.
  int getMolIdx() const { return mol_idx_; }

private:
  /// Flatten ion coordinates into the contiguous layout expected by DeepQMCBridge.
  static std::vector<RealType> flattenIonCoords(const ParticleSet& ions);
  /// Append electron coordinates, substituting activeR(active_iat) for proposed moves.
  static void appendElectronCoords(const ParticleSet& electrons,
                                   std::vector<RealType>& electron_coords,
                                   int active_iat = -1);

  /// Evaluate one DeepQMC batch for all walkers in wfc_list and p_list.
  static DeepQMCBridge::BatchResult evaluateBatch(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                                  int active_iat = -1);
  /// Evaluate a single walker by forwarding through the batched path.
  DeepQMCBridge::BatchResult evaluateOne(const ParticleSet& electrons, int active_iat = -1) const;

  /// Ion ParticleSet used to provide nuclear coordinates.
  const ParticleSet& ions_;
  /// Shared inference bridge so clones can reuse the same Python/JAX model instance.
  std::shared_ptr<const DeepQMCBridge> bridge_;
  /// Molecule index passed through to the DeepQMC model.
  int mol_idx_;
  /// True when ratio() or ratioGrad() has cached a proposed log value.
  bool has_proposed_log_value_ = false;
  /// Cached proposed log value for a particle-by-particle move.
  LogValue proposed_log_value_ = LogValue(0);
};

} // namespace qmcplusplus

#endif
