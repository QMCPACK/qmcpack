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

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"

namespace qmcplusplus
{
class DriftModifierUNR : public DriftModifierBase
{
public:
  using RealType = QMCTraits::RealType;
  using PosType  = QMCTraits::PosType;

  void getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>&) const final;

  void getDrift(RealType tau, const GradType& qf, PosType& drift) const final;

  void getDrifts(RealType tau,
                 const std::vector<ComplexType>& qf,
                 std::vector<ParticleSet::Scalar_t>& drift) const final;

  void getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const final;

  bool parseXML(xmlNodePtr cur) final;

  DriftModifierUNR(RealType a = 1.0) : a_(a) {}

private:
  /// JCP1993 Umrigar et eq. (35) "a" parameter is set to 1.0
  RealType a_;
};

} // namespace qmcplusplus

#endif
