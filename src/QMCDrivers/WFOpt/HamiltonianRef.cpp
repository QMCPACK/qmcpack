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

  void HamiltonianRef::addOperator(OperatorBase& op) { Hrefs_.emplace_back(op); }

  FullPrecRealType HamiltonianRef::evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               std::vector<ValueType>& dlogpsi,
                                               std::vector<ValueType>& dhpsioverpsi,
                                               bool compute_deriv)
  {
    FullPrecRealType LocalEnergy = Hrefs_[0].get().evaluate(P);
    if (compute_deriv)
      for (int i = 1; i < Hrefs_.size(); ++i)
        LocalEnergy += Hrefs_[i].get().evaluateValueAndDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
    else
      for (int i = 1; i < Hrefs_.size(); ++i)
        LocalEnergy += Hrefs_[i].get().evaluate(P);
    return LocalEnergy;
  }

  /// the same evaluate as QMCHamiltonian
  FullPrecRealType HamiltonianRef::evaluate(ParticleSet& P)
  {
    FullPrecRealType LocalEnergy = 0.0;
    for (int i = 0; i < Hrefs_.size(); ++i)
    {
      const auto LocalEnergyComponent = Hrefs_[i].get().evaluate(P);
      if (std::isnan(LocalEnergyComponent))
        APP_ABORT("QMCHamiltonian::evaluate component " + Hrefs_[i].get().myName + " returns NaN\n");
      LocalEnergy += LocalEnergyComponent;
    }
    return LocalEnergy;
  }

} // namespace qmcplusplus
