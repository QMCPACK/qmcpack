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
  if (sknot_ < 2)
    APP_ABORT("Spin knots must be greater than 2\n");

  if (sknot_ % 2 != 0)
    APP_ABORT("Spin knots uses Simpson's rule. Must have even number of knots");

  const size_t num_vars = optvars.num_active_vars;
  dratio_.resize(total_knots_, num_vars);
  dlogpsi_vp_.resize(dlogpsi.size());
  wvec_.resize(total_knots_);

  RealType sold = W.spins[iel];
  buildTotalQuadrature(r, dr, sold);

  //Now we have all the spin and spatial quadrature points acculated to use in evaluation
  //Now we need to obtain dlogpsi and dlogpsi_vp
  if (VP_)
  {
    // Compute ratios with VP
    VP_->makeMovesWithSpin(W, iel, deltaV_, deltaS_, true, iat);
    Psi.evaluateDerivRatios(*VP_, optvars, psiratio_, dratio_);
  }
  else
    for (int iq = 0; iq < total_knots_; iq++)
    {
      PosType posold = W.R[iel];
      W.makeMoveWithSpin(iel, deltaV_[iq], deltaS_[iq]);
      psiratio_[iq] = Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);

      std::fill(dlogpsi_vp_.begin(), dlogpsi_vp_.end(), 0.0);
      Psi.evaluateDerivativesWF(W, optvars, dlogpsi_vp_);
      for (int v = 0; v < dlogpsi_vp_.size(); ++v)
        dratio_(iq, v) = dlogpsi_vp_[v] - dlogpsi[v];

      W.makeMoveWithSpin(iel, -deltaV_[iq], -deltaS_[iq]);
      Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);
    }

  ComplexType pairpot;
  for (int iq = 0; iq < total_knots_; iq++)
  {
    ComplexType lsum;
    for (int il = 0; il < nchannel_; il++)
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
            ldots += lmMatrixElements(l, m1, m2, id) * sMatrixElements(W.spins[iel], W.spins[iel] + deltaS_[iq], id);
          ComplexType cY = std::conj(sphericalHarmonic(l, m2, rrotsgrid_m_[iq % nknot_]));
          msums += Y * cY * ldots;
        }
      }
      lsum += sopp_m_[il]->splint(r) * msums;
    }
    wvec_[iq] = lsum * psiratio_[iq] * spin_quad_weights_[iq];
    pairpot += wvec_[iq];
  }

  BLAS::gemv('N', num_vars, total_knots_, 1.0, dratio_.data(), num_vars, wvec_.data(), 1, 1.0, dhpsioverpsi.data(), 1);

  return std::real(pairpot);
}

} // namespace qmcplusplus
