//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "DistanceTable.h"
#include "CPU/BLAS.hpp"
#include "Utilities/Timer.h"

namespace qmcplusplus
{
NonLocalECPotential::Return_t NonLocalECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                               const opt_variables_type& optvars,
                                                                               const Vector<ValueType>& dlogpsi,
                                                                               Vector<ValueType>& dhpsioverpsi)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG));

  /* evaluating TWF ratio values requires calling prepareGroup
   * In evaluate() we first loop over species and call prepareGroup before looping over all the electrons of a species
   * Here it is not necessary because TWF::evaluateLog has been called and precomputed data is up-to-date
   */
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

/** evaluate the non-local potential of the iat-th ionic center
   * @param W electron configuration
   * @param iat ionic index
   * @param psi trial wavefunction
   * @param optvars optimizables 
   * @param dlogpsi derivatives of the wavefunction at W.R 
   * @param hdpsioverpsi derivatives of Vpp 
   * @param return the non-local component
   *
   * This is a temporary solution which uses TrialWaveFunction::evaluateDerivatives
   * assuming that the distance tables are fully updated for each ratio computation.
   */
NonLocalECPComponent::RealType NonLocalECPComponent::evaluateValueAndDerivatives(ParticleSet& W,
                                                                                 int iat,
                                                                                 TrialWaveFunction& psi,
                                                                                 int iel,
                                                                                 RealType r,
                                                                                 const PosType& dr,
                                                                                 const opt_variables_type& optvars,
                                                                                 const Vector<ValueType>& dlogpsi,
                                                                                 Vector<ValueType>& dhpsioverpsi)
{
  const size_t num_vars = optvars.num_active_vars;
  dratio.resize(nknot, num_vars);
  dlogpsi_vp.resize(dlogpsi.size());

  deltaV.resize(nknot);

  //displacements wrt W.R[iel]
  for (int j = 0; j < nknot; j++)
    deltaV[j] = r * rrotsgrid_m[j] - dr;

  if (VP)
  {
    // Compute ratios with VP
    VP->makeMoves(W, iel, deltaV, true, iat);
    psi.evaluateDerivRatios(*VP, optvars, psiratio, dratio);
  }
  else
  {
    for (int j = 0; j < nknot; j++)
    {
      PosType pos_now = W.R[iel];
      W.makeMove(iel, deltaV[j]);
      psiratio[j] = psi.calcRatio(W, iel);
      psi.acceptMove(W, iel);
      W.acceptMove(iel);

      //use existing methods
      std::fill(dlogpsi_vp.begin(), dlogpsi_vp.end(), 0.0);
      psi.evaluateDerivativesWF(W, optvars, dlogpsi_vp);
      for (int v = 0; v < dlogpsi_vp.size(); ++v)
        dratio(j, v) = dlogpsi_vp[v] - dlogpsi[v];

      W.makeMove(iel, -deltaV[j]);
      psi.calcRatio(W, iel);
      psi.acceptMove(W, iel);
      W.acceptMove(iel);
    }
  }

  for (int j = 0; j < nknot; ++j)
    psiratio[j] *= sgridweight_m[j];

  for (int ip = 0; ip < nchannel; ip++)
    vrad[ip] = nlpp_m[ip]->splint(r) * wgt_angpp_m[ip];

  RealType pairpot(0);
  const RealType rinv = RealType(1) / r;
  // Compute spherical harmonics on grid
  for (int j = 0, jl = 0; j < nknot; j++)
  {
    RealType zz = dot(dr, rrotsgrid_m[j]) * rinv;
    // Forming the Legendre polynomials
    lpol[0]           = 1.0;
    RealType lpolprev = 0.0;
    for (int l = 0; l < lmax; l++)
    {
      lpol[l + 1] = (Lfactor1[l] * zz * lpol[l] - l * lpolprev) * Lfactor2[l];
      lpolprev    = lpol[l];
    }

    RealType lsum = 0.0;
    for (int l = 0; l < nchannel; l++)
      lsum += vrad[l] * lpol[angpp_m[l]];

    wvec[j] = lsum * psiratio[j];
    pairpot += std::real(wvec[j]);
  }

  BLAS::gemv('N', num_vars, nknot, 1.0, dratio.data(), num_vars, wvec.data(), 1, 1.0, dhpsioverpsi.data(), 1);

  return pairpot;
}

} // namespace qmcplusplus
