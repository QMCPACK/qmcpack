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
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{

SOECPotential::Return_t SOECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& optvars,
                                                                   const Vector<ValueType>& dlogpsi,
                                                                   Vector<ValueType>& dhpsioverpsi)
{
  return 0.0;
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
