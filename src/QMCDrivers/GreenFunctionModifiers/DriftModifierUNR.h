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


#ifndef QMCPLUSPLUS_DRIFTMODIFIER_UNR_H
#define QMCPLUSPLUS_DRIFTMODIFIER_UNR_H

/** @file
 *  @brief Declaration of DriftModifierUNR implementing Umrigar-Nightingale-Runge drift scaling.
 */

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"

namespace qmcplusplus
{
/** Implements the Umrigar-Nightingale-Runge (UNR) drift modifier.
 *
 *  Scales the drift velocity according to Eq. (35) of C. J. Umrigar, M. P. Nightingale,
 *  and K. J. Runge, J. Chem. Phys. 99, 2865 (1993):
 *  \f[
 *    \mathbf{d} = \frac{-1 + \sqrt{1 + 2 a \tau v^2}}{a v^2} \mathbf{v}
 *  \f]
 *  where \f$\mathbf{v} = \mathrm{Re}(\nabla \Psi_T / \Psi_T)\f$, \f$v^2 = |\mathbf{v}|^2\f$,
 *  \f$\tau\f$ is the effective timestep, and \f$a\f$ is a tunable positive parameter (default 1.0).
 */
class DriftModifierUNR : public DriftModifierBase
{
public:
  using RealType = QMCTraits::RealType;
  using PosType  = QMCTraits::PosType;

  void getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>& drift) const final;

  void getDrift(RealType tau, const GradType& qf, PosType& drift) const final;

  void getDrifts(RealType tau,
                 const std::vector<ComplexType>& qf,
                 std::vector<ParticleSet::Scalar_t>& drift) const final;

  void getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const final;

  /** parse the UNR parameter "drift_UNR_a" from an XML node
   * @param cur XML node pointer to the QMC driver section
   * @return true if parsing succeeded
   */
  bool parseXML(xmlNodePtr cur) final;

  /** constructor
   * @param a UNR scaling parameter (default 1.0)
   */
  DriftModifierUNR(RealType a = 1.0) : a_(a) {}

private:
  /// JCP1993 Umrigar et eq. (35) "a" parameter is set to 1.0 by default
  RealType a_;
};

} // namespace qmcplusplus

#endif
