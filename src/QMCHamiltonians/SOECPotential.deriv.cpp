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
                                                                     TrialWaveFunction& Psi,
                                                                     int iel,
                                                                     RealType r,
                                                                     const PosType& dr,
                                                                     const opt_variables_type& optvars,
                                                                     const Vector<ValueType>& dlogpsi,
                                                                     Vector<ValueType>& dhpsioverpsi)
{
  if (sknot < 2)
    APP_ABORT("Spin knots must be greater than 2\n");

  if (sknot % 2 != 0)
    APP_ABORT("Spin knots uses Simpson's rule. Must have even number of knots");

  for (int ip = 0; ip < nchannel; ip++)
  {
    vrad[ip] = sopp_m[ip]->splint(r);
  }

  RealType smin(0.0);
  RealType smax(TWOPI);
  RealType dS = (smax - smin) / sknot; //step size for spin

  RealType sold = W.spins[iel];
  ComplexType sint(0.0);

  for (int is = 1; is <= sknot - 1; is += 2)
  {
    RealType snew      = smin + is * dS;
    ComplexType angint = getAngularIntegral(sold, snew, W, Psi, iel, r, dr, iat);
    sint += RealType(4. / 3.) * dS * angint;
  }
  for (int is = 2; is <= sknot - 2; is += 2)
  {
    RealType snew      = smin + is * dS;
    ComplexType angint = getAngularIntegral(sold, snew, W, Psi, iel, r, dr, iat);
    sint += RealType(2. / 3.) * dS * angint;
  }
  sint += RealType(1. / 3.) * dS * getAngularIntegral(sold, smin, W, Psi, iel, r, dr, iat);
  sint += RealType(1. / 3.) * dS * getAngularIntegral(sold, smax, W, Psi, iel, r, dr, iat);

  RealType pairpot = std::real(sint) / TWOPI;
  return pairpot;
}

} // namespace qmcplusplus
