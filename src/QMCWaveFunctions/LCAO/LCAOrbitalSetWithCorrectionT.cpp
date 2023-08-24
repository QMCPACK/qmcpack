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


#include "LCAOrbitalSetWithCorrectionT.h"

namespace qmcplusplus
{
template<typename T>
LCAOrbitalSetWithCorrectionT<T>::LCAOrbitalSetWithCorrectionT(const std::string& my_name,
                                                              ParticleSet& ions,
                                                              ParticleSet& els,
                                                              std::unique_ptr<basis_type>&& bs)
    : SPOSetT<T>(my_name), lcao(my_name + "_modified", std::move(bs)), cusp(ions, els)
{}

template<typename T>
void LCAOrbitalSetWithCorrectionT<T>::setOrbitalSetSize(int norbs)
{
  assert(lcao.getOrbitalSetSize() == norbs && "norbs doesn't agree with lcao!");
  this->OrbitalSetSize = norbs;
  cusp.setOrbitalSetSize(norbs);
}

template<typename T>
std::unique_ptr<SPOSetT<T>> LCAOrbitalSetWithCorrectionT<T>::makeClone() const
{
  return std::make_unique<LCAOrbitalSetWithCorrectionT<T>>(*this);
}

template<typename T>
void LCAOrbitalSetWithCorrectionT<T>::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  lcao.evaluateValue(P, iat, psi);
  cusp.addV(P, iat, psi);
}

template<typename T>
void LCAOrbitalSetWithCorrectionT<T>::evaluateVGL(const ParticleSet& P,
                                                  int iat,
                                                  ValueVector& psi,
                                                  GradVector& dpsi,
                                                  ValueVector& d2psi)
{
  lcao.evaluateVGL(P, iat, psi, dpsi, d2psi);
  cusp.add_vector_vgl(P, iat, psi, dpsi, d2psi);
}

template<typename T>
void LCAOrbitalSetWithCorrectionT<T>::evaluate_notranspose(const ParticleSet& P,
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

template class LCAOrbitalSetWithCorrectionT<double>;
template class LCAOrbitalSetWithCorrectionT<float>;

} // namespace qmcplusplus
