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

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
/// this class implements drift modification
class DriftModifierBase
{
public:
  using RealType    = QMCTraits::RealType;
  using PosType     = QMCTraits::PosType;
  using GradType    = QMCTraits::GradType;
  using ComplexType = QMCTraits::ComplexType;

  /** evaluate a drift with a real force
   * @param tau timestep
   * @param qf quantum force
   * @param drift output
   */
  virtual void getDrift(RealType tau, const GradType& qf, PosType& drift) const = 0;

  virtual void getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const = 0;

  virtual void getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>&) const = 0;

  virtual void getDrifts(RealType tau,
                         const std::vector<ComplexType>& qf,
                         std::vector<ParticleSet::Scalar_t>& drift) const = 0;

  virtual bool parseXML(xmlNodePtr cur) { return true; }

  virtual ~DriftModifierBase() {}

protected:
  // modifer name
  std::string ClassName;
};

} // namespace qmcplusplus

#endif
