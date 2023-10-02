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
LCAOrbitalSetWithCorrection::LCAOrbitalSetWithCorrection(const std::string& my_name,
                                                         ParticleSet& ions,
                                                         ParticleSet& els,
                                                         std::unique_ptr<basis_type>&& bs)
    : SPOSet(my_name), lcao(my_name + "_modified", std::move(bs)), cusp(ions, els)
{}

void LCAOrbitalSetWithCorrection::setOrbitalSetSize(int norbs)
{
  assert(lcao.getOrbitalSetSize() == norbs && "norbs doesn't agree with lcao!");
  OrbitalSetSize = norbs;
  cusp.setOrbitalSetSize(norbs);
}


std::unique_ptr<SPOSet> LCAOrbitalSetWithCorrection::makeClone() const
{
  return std::make_unique<LCAOrbitalSetWithCorrection>(*this);
}

void LCAOrbitalSetWithCorrection::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  lcao.evaluateValue(P, iat, psi);
  cusp.addV(P, iat, psi);
}

void LCAOrbitalSetWithCorrection::evaluateVGL(const ParticleSet& P,
                                              int iat,
                                              ValueVector& psi,
                                              GradVector& dpsi,
                                              ValueVector& d2psi)
{
  lcao.evaluateVGL(P, iat, psi, dpsi, d2psi);
  cusp.add_vector_vgl(P, iat, psi, dpsi, d2psi);
}

void LCAOrbitalSetWithCorrection::evaluate_notranspose(const ParticleSet& P,
                                                       int first,
                                                       int last,
                                                       ValueMatrix& logdet,
                                                       GradMatrix& dlogdet,
                                                       ValueMatrix& d2logdet)
{
  lcao.evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
  for (size_t i = 0, iat = first; iat < last; i++, iat++)
    cusp.add_vgl(P, iat, i, logdet, dlogdet, d2logdet);
}

} // namespace qmcplusplus
