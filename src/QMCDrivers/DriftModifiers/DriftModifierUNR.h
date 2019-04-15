//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DRIFTMODIFIER_UNR_H
#define QMCPLUSPLUS_DRIFTMODIFIER_UNR_H

#include "QMCDrivers/DriftModifiers/DriftModifierBase.h"

namespace qmcplusplus
{
class DriftModifierUNR : public DriftModifierBase
{
public:
  using RealType = QMCTraits::RealType;
  using PosType  = QMCTraits::PosType;

  void getScaledDrift(const ParticleSet& P,
                      const TrialWaveFunction& W,
                      const QMCHamiltonian& H,
                      RealType tau,
                      const GradType& qf,
                      int iat,
                      PosType& drift) const final;

  bool parseXML(xmlNodePtr cur) final { return true; };
};

} // namespace qmcplusplus

#endif
