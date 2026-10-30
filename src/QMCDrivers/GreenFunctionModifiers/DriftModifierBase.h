//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DRIFTMODIFIER_BASE_H
#define QMCPLUSPLUS_DRIFTMODIFIER_BASE_H

/** @file
 *  @brief Declaration of the base class DriftModifierBase for drift velocity modifiers.
 */

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TWFGrads.hpp"
#include "TauParams.hpp"

namespace qmcplusplus
{
/** Base class for drift modification schemes in QMC drivers.
 *
 *  Drift modifiers scale or transform the bare wavefunction gradient (quantum force)
 *  to compute the displacement drift vector used in Langevin/diffusion moves, preventing
 *  divergences near nodal surfaces where gradients become large.
 */
class DriftModifierBase
{
public:
  using RealType    = QMCTraits::RealType;
  using PosType     = QMCTraits::PosType;
  using GradType    = QMCTraits::GradType;
  using ComplexType = QMCTraits::ComplexType;

  /** evaluate a spatial drift for a single particle from a wavefunction gradient
   * @param tau effective timestep (typically tau / mass)
   * @param qf spatial wavefunction gradient (quantum force)
   * @param[out] drift output scaled spatial drift vector
   */
  virtual void getDrift(RealType tau, const GradType& qf, PosType& drift) const = 0;

  /** evaluate a scalar (spin) drift for a single particle from a spin gradient
   * @param tau effective spin timestep (typically tau / (mass * spin_mass))
   * @param qf complex spin gradient
   * @param[out] drift output scaled scalar spin drift
   */
  virtual void getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const = 0;

  /** evaluate spatial drifts for a vector of particle gradients
   * @param tau effective timestep (typically tau / mass)
   * @param qf vector of spatial wavefunction gradients
   * @param[out] drift output vector of scaled spatial drift vectors
   */
  virtual void getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>& drift) const = 0;

  /** evaluate scalar (spin) drifts for a vector of spin gradients
   * @param tau effective spin timestep (typically tau / (mass * spin_mass))
   * @param qf vector of complex spin gradients
   * @param[out] drift output vector of scaled scalar spin drifts
   */
  virtual void getDrifts(RealType tau,
                         const std::vector<ComplexType>& qf,
                         std::vector<ParticleSet::Scalar_t>& drift) const = 0;

  /** evaluate drifts for all coordinates (spatial and optionally spin) given TauParams and TWFGrads
   * @tparam CT coordinate type (CoordsType::POS or CoordsType::POS_SPIN)
   * @param taus effective timestep parameters for spatial and spin coordinates
   * @param qf trial wavefunction gradients for spatial and spin coordinates
   * @param[out] drifts output drift coordinates
   */
  template<CoordsType CT>
  void getDrifts(const TauParams<RealType, CT>& taus, const TWFGrads<CT>& qf, MCCoords<CT>& drifts) const;

  /** parse drift modifier parameters from an XML node
   * @param cur XML node pointer to the QMC driver section
   * @return true if parsing succeeded
   */
  virtual bool parseXML(xmlNodePtr cur) { return true; }

  virtual ~DriftModifierBase() {}

protected:
  /// modifier name
  std::string ClassName;
};

template<CoordsType CT>
void DriftModifierBase::getDrifts(const TauParams<RealType, CT>& taus, const TWFGrads<CT>& qf, MCCoords<CT>& drifts) const
{
  getDrifts(taus.tauovermass, qf.grads_positions, drifts.positions);
  if constexpr (CT == CoordsType::POS_SPIN)
    getDrifts(taus.spin_tauovermass, qf.grads_spins, drifts.spins);
}

} // namespace qmcplusplus

#endif
