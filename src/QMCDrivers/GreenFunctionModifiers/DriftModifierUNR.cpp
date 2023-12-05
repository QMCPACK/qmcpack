//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "DriftModifierUNR.h"
#include "OhmmsData/ParameterSet.h"
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
void DriftModifierUNR::getDrift(RealType tau, const GradType& qf, PosType& drift) const
{
  // convert the complex WF gradient to real
  convertToReal(qf, drift);
  RealType vsq = dot(drift, drift);
  RealType sc  = vsq < std::numeric_limits<RealType>::epsilon()
      ? tau
      : ((-1.0 + std::sqrt(1.0 + 2.0 * a_ * tau * vsq)) / (a_ * vsq));
  //Apply the umrigar scaling to drift.
  drift *= sc;
}

void DriftModifierUNR::getDrift(RealType tau, const ComplexType& qf, ParticleSet::Scalar_t& drift) const
{
  // convert the complex WF gradient to real
  convertToReal(qf, drift);
  RealType vsq = drift * drift;
  RealType sc  = vsq < std::numeric_limits<RealType>::epsilon()
      ? tau
      : ((-1.0 + std::sqrt(1.0 + 2.0 * a_ * tau * vsq)) / (a_ * vsq));
  //Apply the umrigar scaling to drift.
  drift *= sc;
}

void DriftModifierUNR::getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>& drift) const
{
  for (int i = 0; i < qf.size(); ++i)
    getDrift(tau, qf[i], drift[i]);
}

void DriftModifierUNR::getDrifts(RealType tau,
                                 const std::vector<ComplexType>& qf,
                                 std::vector<ParticleSet::Scalar_t>& drift) const
{
  for (int i = 0; i < qf.size(); ++i)
    getDrift(tau, qf[i], drift[i]);
}

bool DriftModifierUNR::parseXML(xmlNodePtr cur)
{
  ParameterSet m_param;
  m_param.add(a_, "drift_UNR_a");
  m_param.put(cur);
  app_log() << "  Set drift_modifier UNR parameter a = " << a_ << std::endl;
  return true;
}

} // namespace qmcplusplus
