//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Jeongnim Kim and QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "kSpaceJastrow.h"

#include "CPU/e2iphi.h" // for eval_e2iphi

namespace qmcplusplus
{

template<class T>
T kSpaceJastrow::do_ratioT(ParticleSet& P, int iat)
{
  RealType J1new(0.0), J1old(0.0), J2new(0.0), J2old(0.0);
  const PosType &rnew(P.getActivePos()), &rold(P.R[iat]);
  // Compute one-body contribution
  int nOne = OneBodyGvecs.size();
  for (int i = 0; i < nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rnew);
  eval_e2iphi(OneBodyPhase, OneBody_e2iGr);
  for (int i = 0; i < nOne; i++)
    J1new += Prefactor * real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  for (int i = 0; i < nOne; i++)
    OneBodyPhase[i] = dot(OneBodyGvecs[i], rold);
  eval_e2iphi(OneBodyPhase, OneBody_e2iGr);
  for (int i = 0; i < nOne; i++)
    J1old += Prefactor * real(OneBodyCoefs[i] * qmcplusplus::conj(OneBody_e2iGr[i]));
  // Now, do two-body part
  int nTwo = TwoBodyGvecs.size();
  for (int i = 0; i < nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rold);
  eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_old);
  for (int i = 0; i < nTwo; i++)
    TwoBodyPhase[i] = dot(TwoBodyGvecs[i], rnew);
  eval_e2iphi(TwoBodyPhase, TwoBody_e2iGr_new);
  for (int i = 0; i < nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i];
    J2old += Prefactor * TwoBodyCoefs[i] * std::norm(rho_G);
  }
  for (int i = 0; i < nTwo; i++)
  {
    ComplexType rho_G = TwoBody_rhoG[i] + TwoBody_e2iGr_new[i] - TwoBody_e2iGr_old[i];
    J2new += Prefactor * TwoBodyCoefs[i] * std::norm(rho_G);
  }
  return std::exp(static_cast<T>(J1new + J2new - (J1old + J2old)));
}

} // namespace qmcplusplus