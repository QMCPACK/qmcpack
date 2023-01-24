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
#include "Numerics/Ylm.h"
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

  int total_knots = nknot * (sknot + 1);
  deltaV.resize(total_knots);
  deltaS.resize(total_knots);
  wvec.resize(total_knots);
  spin_quad_weights.resize(total_knots);
  psiratio.resize(total_knots);
  const size_t num_vars = optvars.num_active_vars;
  dratio.resize(total_knots, num_vars);
  dlogpsi_vp.resize(dlogpsi.size());

  auto addAngularQuadraturePoints = [&](const int is, const RealType r, const PosType& dr, RealType ds, RealType swt) {
    for (int iq = 0; iq < nknot; iq++)
    {
      int offset     = is * nknot + iq;
      deltaV[offset] = r * rrotsgrid_m[iq] - dr;
      deltaS[offset] = ds;
      //spin_norm is 1/(2pi), angular_norm is 4pi, so total norm is 2
      spin_quad_weights[offset] = 2.0 * swt * sgridweight_m[iq];
    }
  };

  RealType smin(0.0);
  RealType smax(TWOPI);
  RealType dS = (smax - smin) / sknot; //step size for spin

  RealType sold = W.spins[iel];
  ComplexType pairpot;
  for (int is = 1; is <= sknot - 1; is += 2)
  {
    RealType snew = smin + is * dS;
    addAngularQuadraturePoints(is, r, dr, snew - sold, RealType(4. / 3.) * dS);
  }
  for (int is = 2; is <= sknot - 2; is += 2)
  {
    RealType snew = smin + is * dS;
    addAngularQuadraturePoints(is, r, dr, snew - sold, RealType(2. / 3.) * dS);
  }
  addAngularQuadraturePoints(0, r, dr, smin - sold, RealType(1. / 3.) * dS);
  addAngularQuadraturePoints(sknot, r, dr, smax - sold, RealType(1. / 3.) * dS);

  //Now we have all the spin and spatial quadrature points acculated to use in evaluation
  //Now we need to obtain dlogpsi and dlogpsi_vp
  if (VP)
  {
    // Compute ratios with VP
    VP->makeMovesWithSpin(W, iel, deltaV, deltaS, true, iat);
    Psi.evaluateDerivRatios(*VP, optvars, psiratio, dratio);
  }
  else
  {
    for (int iq = 0; iq < total_knots; iq++)
    {
      PosType pos_now   = W.R[iel];
      RealType spin_now = W.spins[iel];
      W.makeMoveWithSpin(iel, deltaV[iq], deltaS[iq]);
      psiratio[iq] = Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);

      //use existing methods
      std::fill(dlogpsi_vp.begin(), dlogpsi_vp.end(), 0.0);
      Psi.evaluateDerivativesWF(W, optvars, dlogpsi_vp);
      for (int v = 0; v < dlogpsi_vp.size(); ++v)
        dratio(iq, v) = dlogpsi_vp[v] - dlogpsi[v];

      W.makeMoveWithSpin(iel, -deltaV[iq], -deltaS[iq]);
      Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);
    }
  }

  for (int iq = 0; iq < total_knots; iq++)
  {
    ComplexType lsum;
    for (int il = 0; il < nchannel; il++)
    {
      int l = il + 1;
      ComplexType msums;
      for (int m1 = -l; m1 <= l; m1++)
      {
        ComplexType Y = sphericalHarmonic(l, m1, dr);
        for (int m2 = -l; m2 <= l; m2++)
        {
          ComplexType ldots;
          for (int id = 0; id < 3; id++)
            ldots += lmMatrixElements(l, m1, m2, id) * sMatrixElements(W.spins[iel], W.spins[iel] + deltaS[iq], id);
          ComplexType cY = std::conj(sphericalHarmonic(l, m2, rrotsgrid_m[iq % (sknot + 1)]));
          msums += Y * cY * ldots;
        }
      }
      lsum += sopp_m[il]->splint(r) * msums;
    }
    pairpot += psiratio[iq] * spin_quad_weights[iq];
  }

  return std::real(pairpot);
}

} // namespace qmcplusplus
