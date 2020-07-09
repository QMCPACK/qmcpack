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


#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
void DriftModifierUNR::getDrift(RealType tau, const GradType& qf, PosType& drift) const
{
  // convert the complex WF gradient to real
  convert(qf, drift);
#ifndef NDEBUG
  PosType debug_drift = drift;
#endif
  RealType vsq = dot(drift, drift);
  RealType sc  = (std::isnan(vsq) || vsq < std::numeric_limits<RealType>::epsilon())
      ? tau
      : ((-1.0 + std::sqrt(1.0 + 2.0 * a_ * tau * vsq)) / (a_ * vsq));
  //Apply the umrigar scaling to drift.
  drift *= sc;
#ifndef NDEBUG
  for(int i = 0; i < drift.size(); ++i)
    assert(!std::isnan(drift[i]));
#endif
}

void DriftModifierUNR::getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const
{
  // convert the complex WF gradient to real
  convert(qf, drift);
  RealType vsq = drift * drift;
  RealType sc  = (std::isnan(vsq) || vsq < std::numeric_limits<RealType>::epsilon())
      ? tau
      : ((-1.0 + std::sqrt(1.0 + 2.0 * a_ * tau * vsq)) / (a_ * vsq));
  //Apply the umrigar scaling to drift.
  drift *= sc;
}

void DriftModifierUNR::getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>& drift) const
{
  for (int i = 0; i < qf.size(); ++i)
  {
    getDrift(tau, qf[i], drift[i]);
  }
}

bool DriftModifierUNR::parseXML(xmlNodePtr cur)
{
  ParameterSet m_param;
  m_param.add(a_, "drift_UNR_a", "double");
  m_param.put(cur);
  app_log() << "  Set drift_modifier UNR parameter a = " << a_ << std::endl;
  return true;
}

} // namespace qmcplusplus
