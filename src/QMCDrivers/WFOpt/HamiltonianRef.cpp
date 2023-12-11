//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "HamiltonianRef.h"

namespace qmcplusplus
{
using FullPrecRealType = HamiltonianRef::FullPrecRealType;

HamiltonianRef::HamiltonianRef(const RefVector<OperatorBase> refs) : Hrefs_(refs) {}

FullPrecRealType HamiltonianRef::evaluateValueAndDerivatives(ParticleSet& P,
                                                             const opt_variables_type& optvars,
                                                             Vector<ValueType>& dlogpsi,
                                                             Vector<ValueType>& dhpsioverpsi)
{
  FullPrecRealType LocalEnergy(0);
  for (OperatorBase& Href : Hrefs_)
    LocalEnergy += Href.evaluateValueAndDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  return LocalEnergy;
}

/// the same evaluate as QMCHamiltonian
FullPrecRealType HamiltonianRef::evaluate(ParticleSet& P)
{
  FullPrecRealType LocalEnergy = 0.0;
  for (int i = 0; i < Hrefs_.size(); ++i)
  {
    const auto LocalEnergyComponent = Hrefs_[i].get().evaluate(P);
    if (qmcplusplus::isnan(LocalEnergyComponent))
      APP_ABORT("HamiltonianRef::evaluate component " + Hrefs_[i].get().getName() + " returns NaN\n");
    LocalEnergy += LocalEnergyComponent;
  }
  return LocalEnergy;
}

} // namespace qmcplusplus
