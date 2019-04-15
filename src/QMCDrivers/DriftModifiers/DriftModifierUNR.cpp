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


#include "QMCDrivers/DriftModifiers/DriftModifierUNR.h"

namespace qmcplusplus
{
void DriftModifierUNR::getScaledDrift(const ParticleSet& P,
                                      const TrialWaveFunction& W,
                                      const QMCHamiltonian& H,
                                      RealType tau,
                                      GradType qf,
                                      int iat,
                                      PosType drift) const
{
  // convert the complex WF gradient to real and temporarily store in drift.
  convert(qf, drift);
  RealType vsq = dot(drift, drift);
  RealType sc =
      (vsq < std::numeric_limits<RealType>::epsilon()) ? tau : ((-1.0 + std::sqrt(1.0 + 2.0 * tau * vsq)) / vsq);
  //Apply the umrigar scaled drift.
  drift *= sc;
}

} // namespace qmcplusplus
