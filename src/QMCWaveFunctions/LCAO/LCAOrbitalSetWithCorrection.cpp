//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include "LCAOrbitalSetWithCorrection.h"

namespace qmcplusplus
{
LCAOrbitalSetWithCorrection::LCAOrbitalSetWithCorrection(ParticleSet& ions,
                                                         ParticleSet& els,
                                                         std::unique_ptr<basis_type>&& bs,
                                                         bool optimize)
    : LCAOrbitalSet(std::move(bs), optimize), cusp(ions, els)
{}

void LCAOrbitalSetWithCorrection::setOrbitalSetSize(int norbs)
{
  LCAOrbitalSet::setOrbitalSetSize(norbs);
  cusp.setBasisSetSize(norbs);
}


std::unique_ptr<SPOSet> LCAOrbitalSetWithCorrection::makeClone() const
{
  return std::make_unique<LCAOrbitalSetWithCorrection>(*this);
}

void LCAOrbitalSetWithCorrection::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  LCAOrbitalSet::evaluateValue(P, iat, psi);
  cusp.addV(P, iat, psi.data());
}

void LCAOrbitalSetWithCorrection::evaluateVGL(const ParticleSet& P,
                                              int iat,
                                              ValueVector& psi,
                                              GradVector& dpsi,
                                              ValueVector& d2psi)
{
  LCAOrbitalSet::evaluateVGL(P, iat, psi, dpsi, d2psi);
  cusp.add_vector_vgl(P, iat, psi, dpsi, d2psi);
}

void LCAOrbitalSetWithCorrection::evaluateVGH(const ParticleSet& P,
                                              int iat,
                                              ValueVector& psi,
                                              GradVector& dpsi,
                                              HessVector& grad_grad_psi)
{
  APP_ABORT("LCAOrbitalSetWithCorrection::evaluate with HessVector not implemented");
}

void LCAOrbitalSetWithCorrection::evaluate_notranspose(const ParticleSet& P,
                                                       int first,
                                                       int last,
                                                       ValueMatrix& logdet,
                                                       GradMatrix& dlogdet,
                                                       ValueMatrix& d2logdet)
{
  LCAOrbitalSet::evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
  for (size_t i = 0, iat = first; iat < last; i++, iat++)
    cusp.add_vgl(P, iat, i, logdet, dlogdet, d2logdet);
}

void LCAOrbitalSetWithCorrection::evaluate_notranspose(const ParticleSet& P,
                                                       int first,
                                                       int last,
                                                       ValueMatrix& logdet,
                                                       GradMatrix& dlogdet,
                                                       HessMatrix& grad_grad_logdet)
{
  APP_ABORT("LCAOrbitalSetWithCorrection::evaluate_notranspose with HessMatrix not implemented");
}

void LCAOrbitalSetWithCorrection::evaluate_notranspose(const ParticleSet& P,
                                                       int first,
                                                       int last,
                                                       ValueMatrix& logdet,
                                                       GradMatrix& dlogdet,
                                                       HessMatrix& grad_grad_logdet,
                                                       GGGMatrix& grad_grad_grad_logdet)
{
  APP_ABORT("LCAOrbitalSetWithCorrection::evaluate_notranspose with GGGMatrix not implemented");
}

void LCAOrbitalSetWithCorrection::evaluateThirdDeriv(const ParticleSet& P,
                                                     int first,
                                                     int last,
                                                     GGGMatrix& grad_grad_grad_logdet)
{
  APP_ABORT("LCAOrbitalSetWithCorrection::evaluateThirdDeriv not implemented");
}
} // namespace qmcplusplus
