//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCHamiltonians/SOECPComponent.h"
#include "QMCHamiltonians/SOECPotential.h"
#include "DistanceTable.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{

SOECPotential::Return_t SOECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& optvars,
                                                                   const Vector<ValueType>& dlogpsi,
                                                                   Vector<ValueType>& dhpsioverpsi)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG));

  const auto& myTable = P.getDistTableAB(myTableIndex);
  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.getDistRow(jel);
    const auto& displ = myTable.getDisplRow(jel);
    for (int iat = 0; iat < NumIons; iat++)
      if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        value_ += PP[iat]->evaluateValueAndDerivatives(P, iat, Psi, jel, dist[iat], -displ[iat], optvars, dlogpsi,
                                                       dhpsioverpsi);
  }
  return value_;
}

SOECPComponent::RealType SOECPComponent::evaluateValueAndDerivatives(ParticleSet& W,
                                                                     int iat,
                                                                     TrialWaveFunction& psi,
                                                                     int iel,
                                                                     RealType r,
                                                                     const PosType& dr,
                                                                     const opt_variables_type& optvars,
                                                                     const Vector<ValueType>& dlogpsi,
                                                                     Vector<ValueType>& dhpsioverpsi)
{
  return 0.0;
}

} // namespace qmcplusplus
