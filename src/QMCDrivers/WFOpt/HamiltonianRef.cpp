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
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif

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
      APP_ABORT("HamiltonianRef::evaluate component " + Hrefs_[i].get().getName() + " returns NaN\n");
    LocalEnergy += LocalEnergyComponent;
  }
  return LocalEnergy;
}

int HamiltonianRef::addObservables(ParticleSet& P)
{
  OperatorBase::PropertySetType Observables;
  //ParticleSet::mcObservables (large data, e.g. density) are accumulated while evaluations
  P.Collectables.clear();
  P.Collectables.rewind();
  for (int i = 0; i < Hrefs_.size(); ++i)
    Hrefs_[i].get().addObservables(Observables, P.Collectables);
  const int myIndex = P.PropertyList.add(Observables.Names[0]);
  for (int i = 1; i < Observables.size(); ++i)
    P.PropertyList.add(Observables.Names[i]);
  P.Collectables.size();
  app_log() << "\n  QMCHamiltonian::add2WalkerProperty added"
            << "\n    " << Observables.size() << " to P::PropertyList "
            << "\n    " << P.Collectables.size() << " to P::Collectables "
            << "\n    starting Index of the observables in P::PropertyList = " << myIndex << std::endl;
  return Observables.size();
}

#ifdef QMC_CUDA
void HamiltonianRef::evaluate(MCWalkerConfiguration& W, std::vector<RealType>& energyVector)
{
  using WP      = WalkerProperties::Indexes;
  auto& walkers = W.WalkerList;
  const int nw  = walkers.size();
  std::vector<FullPrecRealType> LocalEnergyVector(nw, 0.0);
  for (int i = 0; i < Hrefs_.size(); ++i)
    Hrefs_[i].get().addEnergy(W, LocalEnergyVector);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[WP::LOCALENERGY] = LocalEnergyVector[iw];
    walkers[iw]->getPropertyBase()[WP::LOCALPOTENTIAL] =
        LocalEnergyVector[iw] - walkers[iw]->getPropertyBase()[WP::NUMPROPERTIES];
    energyVector[iw] = LocalEnergyVector[iw];
  }
}
#endif

} // namespace qmcplusplus
