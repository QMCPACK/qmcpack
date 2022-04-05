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
#ifndef NDEBUG
  PosType debug_drift = drift;
#endif
  RealType vsq = dot(drift, drift);
  RealType sc  = vsq < std::numeric_limits<RealType>::epsilon()
      ? tau
      : ((-1.0 + std::sqrt(1.0 + 2.0 * a_ * tau * vsq)) / (a_ * vsq));
  //Apply the umrigar scaling to drift.
  drift *= sc;
#ifndef NDEBUG
  // Why is this in NDEBUG: at least for gnu std::isnan may be no-op if NDEBUG is defined.
  // Generally we hope that this would only occur as the result of bad input
  // which would hopefully be the result of development time error and
  // therefore caught when run in a Debug build.
  if (std::isnan(vsq))
  {
    std::ostringstream error_message;
    for (int i = 0; i < drift.size(); ++i)
    {
      if (std::isnan(drift[i]))
      {
        error_message << "drift[" << i << "] is nan, vsq (" << vsq << ") sc (" << sc << ")\n";
        break;
      }
    }
    throw std::runtime_error(error_message.str());
  }
#endif
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
#ifndef NDEBUG
  // Why is this in NDEBUG: at least for gnu std::isnan may be no-op if NDEBUG is defined.
  // Generally we hope that this would only occur as the result of bad input
  // which would hopefully be the result of development time error and
  // therefore caught when run in a Debug build.
  if (std::isnan(vsq))
  {
    std::ostringstream error_message;
    if (std::isnan(drift))
    {
      error_message << "drift is nan, vsq (" << vsq << ") sc (" << sc << ")\n";
    }
    else
    {
      error_message << "vsq is nan but drift is " << drift << ", unexpected, investigate.\n";
    }
    throw std::runtime_error(error_message.str());
  }
#endif
}

void DriftModifierUNR::getDrifts(RealType tau, const std::vector<GradType>& qf, std::vector<PosType>& drift) const
{
  for (int i = 0; i < qf.size(); ++i)
  {
    getDrift(tau, qf[i], drift[i]);
  }
}

void DriftModifierUNR::getDrifts(RealType tau,
                                 const std::vector<ComplexType>& qf,
                                 std::vector<ParticleSet::Scalar_t>& drift) const
{
  for (int i = 0; i < qf.size(); ++i)
  {
    getDrift(tau, qf[i], drift[i]);
  }
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
