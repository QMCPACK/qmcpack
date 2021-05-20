//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:  Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "SpinorSet.h"

namespace qmcplusplus
{
SpinorSet::SpinorSet() : SPOSet(), className("SpinorSet"), spo_up(nullptr), spo_dn(nullptr) {}

void SpinorSet::set_spos(std::unique_ptr<SPOSet>&& up, std::unique_ptr<SPOSet>&& dn)
{
  //Sanity check for input SPO's.  They need to be the same size or
  IndexType spo_size_up   = up->getOrbitalSetSize();
  IndexType spo_size_down = dn->getOrbitalSetSize();

  if (spo_size_up != spo_size_down)
    APP_ABORT("SpinorSet::set_spos(...):  up and down SPO components have different sizes.");

  setOrbitalSetSize(spo_size_up);

  spo_up = std::move(up);
  spo_dn = std::move(dn);

  psi_work_up.resize(OrbitalSetSize);
  psi_work_down.resize(OrbitalSetSize);

  dpsi_work_up.resize(OrbitalSetSize);
  dpsi_work_down.resize(OrbitalSetSize);

  d2psi_work_up.resize(OrbitalSetSize);
  d2psi_work_down.resize(OrbitalSetSize);
}

int SpinorSet::getBasisSetSize() const
{
  IndexType basis_up = spo_up->getBasisSetSize();
  IndexType basis_dn = spo_dn->getBasisSetSize();
  if (basis_up != basis_dn)
    APP_ABORT("SpinorSet::getBasisSetSize(): up and down basis size not equal.");
  return basis_up;
}

void SpinorSet::resetParameters(const opt_variables_type& optVariables){};

void SpinorSet::setOrbitalSetSize(int norbs) { OrbitalSetSize = norbs; };


void SpinorSet::evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  psi_work_up   = 0.0;
  psi_work_down = 0.0;

  spo_up->evaluateValue(P, iat, psi_work_up);
  spo_dn->evaluateValue(P, iat, psi_work_down);

  ParticleSet::Scalar_t s = P.activeSpin(iat);

  RealType coss(0.0), sins(0.0);

  coss = std::cos(s);
  sins = std::sin(s);

  //This is only supported in the complex build, so ValueType is some complex number depending on the precision.
  ValueType eis(coss, sins);
  ValueType emis(coss, -sins);

  psi = eis * psi_work_up + emis * psi_work_down;
}

void SpinorSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  psi_work_up     = 0.0;
  psi_work_down   = 0.0;
  dpsi_work_up    = 0.0;
  dpsi_work_down  = 0.0;
  d2psi_work_up   = 0.0;
  d2psi_work_down = 0.0;

  spo_up->evaluateVGL(P, iat, psi_work_up, dpsi_work_up, d2psi_work_up);
  spo_dn->evaluateVGL(P, iat, psi_work_down, dpsi_work_down, d2psi_work_down);

  ParticleSet::Scalar_t s = P.activeSpin(iat);

  RealType coss(0.0), sins(0.0);

  coss = std::cos(s);
  sins = std::sin(s);

  ValueType eis(coss, sins);
  ValueType emis(coss, -sins);

  psi   = eis * psi_work_up + emis * psi_work_down;
  dpsi  = eis * dpsi_work_up + emis * dpsi_work_down;
  d2psi = eis * d2psi_work_up + emis * d2psi_work_down;
}

void SpinorSet::evaluateVGL_spin(const ParticleSet& P,
                                 int iat,
                                 ValueVector_t& psi,
                                 GradVector_t& dpsi,
                                 ValueVector_t& d2psi,
                                 ValueVector_t& dspin)
{
  psi_work_up     = 0.0;
  psi_work_down   = 0.0;
  dpsi_work_up    = 0.0;
  dpsi_work_down  = 0.0;
  d2psi_work_up   = 0.0;
  d2psi_work_down = 0.0;

  spo_up->evaluateVGL(P, iat, psi_work_up, dpsi_work_up, d2psi_work_up);
  spo_dn->evaluateVGL(P, iat, psi_work_down, dpsi_work_down, d2psi_work_down);

  ParticleSet::Scalar_t s = P.activeSpin(iat);

  RealType coss(0.0), sins(0.0);

  coss = std::cos(s);
  sins = std::sin(s);

  ValueType eis(coss, sins);
  ValueType emis(coss, -sins);
  ValueType eye(0, 1.0);

  psi   = eis * psi_work_up + emis * psi_work_down;
  dpsi  = eis * dpsi_work_up + emis * dpsi_work_down;
  d2psi = eis * d2psi_work_up + emis * d2psi_work_down;
  dspin = eye * (eis * psi_work_up - emis * psi_work_down);
}

void SpinorSet::evaluate_notranspose(const ParticleSet& P,
                                     int first,
                                     int last,
                                     ValueMatrix_t& logdet,
                                     GradMatrix_t& dlogdet,
                                     ValueMatrix_t& d2logdet)
{
  IndexType nelec = P.getTotalNum();

  logpsi_work_up.resize(nelec, OrbitalSetSize);
  logpsi_work_down.resize(nelec, OrbitalSetSize);

  dlogpsi_work_up.resize(nelec, OrbitalSetSize);
  dlogpsi_work_down.resize(nelec, OrbitalSetSize);

  d2logpsi_work_up.resize(nelec, OrbitalSetSize);
  d2logpsi_work_down.resize(nelec, OrbitalSetSize);

  spo_up->evaluate_notranspose(P, first, last, logpsi_work_up, dlogpsi_work_up, d2logpsi_work_up);
  spo_dn->evaluate_notranspose(P, first, last, logpsi_work_down, dlogpsi_work_down, d2logpsi_work_down);


  for (int iat = 0; iat < nelec; iat++)
  {
    ParticleSet::Scalar_t s = P.activeSpin(iat);

    RealType coss(0.0), sins(0.0);

    coss = std::cos(s);
    sins = std::sin(s);

    ValueType eis(coss, sins);
    ValueType emis(coss, -sins);

    for (int no = 0; no < OrbitalSetSize; no++)
    {
      logdet(iat, no)   = eis * logpsi_work_up(iat, no) + emis * logpsi_work_down(iat, no);
      dlogdet(iat, no)  = eis * dlogpsi_work_up(iat, no) + emis * dlogpsi_work_down(iat, no);
      d2logdet(iat, no) = eis * d2logpsi_work_up(iat, no) + emis * d2logpsi_work_down(iat, no);
    }
  }
}

void SpinorSet::evaluate_notranspose_spin(const ParticleSet& P,
                                          int first,
                                          int last,
                                          ValueMatrix_t& logdet,
                                          GradMatrix_t& dlogdet,
                                          ValueMatrix_t& d2logdet,
                                          ValueMatrix_t& dspinlogdet)
{
  IndexType nelec = P.getTotalNum();

  logpsi_work_up.resize(nelec, OrbitalSetSize);
  logpsi_work_down.resize(nelec, OrbitalSetSize);

  dlogpsi_work_up.resize(nelec, OrbitalSetSize);
  dlogpsi_work_down.resize(nelec, OrbitalSetSize);

  d2logpsi_work_up.resize(nelec, OrbitalSetSize);
  d2logpsi_work_down.resize(nelec, OrbitalSetSize);

  spo_up->evaluate_notranspose(P, first, last, logpsi_work_up, dlogpsi_work_up, d2logpsi_work_up);
  spo_dn->evaluate_notranspose(P, first, last, logpsi_work_down, dlogpsi_work_down, d2logpsi_work_down);


  for (int iat = 0; iat < nelec; iat++)
  {
    ParticleSet::Scalar_t s = P.activeSpin(iat);

    RealType coss(0.0), sins(0.0);

    coss = std::cos(s);
    sins = std::sin(s);

    ValueType eis(coss, sins);
    ValueType emis(coss, -sins);
    ValueType eye(0, 1.0);

    for (int no = 0; no < OrbitalSetSize; no++)
    {
      logdet(iat, no)      = eis * logpsi_work_up(iat, no) + emis * logpsi_work_down(iat, no);
      dlogdet(iat, no)     = eis * dlogpsi_work_up(iat, no) + emis * dlogpsi_work_down(iat, no);
      d2logdet(iat, no)    = eis * d2logpsi_work_up(iat, no) + emis * d2logpsi_work_down(iat, no);
      dspinlogdet(iat, no) = eye * (eis * logpsi_work_up(iat, no) - emis * logpsi_work_down(iat, no));
    }
  }
}


void SpinorSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector_t& psi, ValueVector_t& dpsi)
{
  psi_work_up   = 0.0;
  psi_work_down = 0.0;

  spo_up->evaluateValue(P, iat, psi_work_up);
  spo_dn->evaluateValue(P, iat, psi_work_down);

  ParticleSet::Scalar_t s = P.activeSpin(iat);

  RealType coss(0.0), sins(0.0);

  coss = std::cos(s);
  sins = std::sin(s);

  ValueType eis(coss, sins);
  ValueType emis(coss, -sins);
  ValueType eye(0, 1.0);

  psi  = eis * psi_work_up + emis * psi_work_down;
  dpsi = eye * (eis * psi_work_up - emis * psi_work_down);
}

SPOSet* SpinorSet::makeClone() const
{
  SpinorSet* myclone = new SpinorSet();
  std::unique_ptr<SPOSet> cloneup(spo_up->makeClone());
  std::unique_ptr<SPOSet> clonedn(spo_dn->makeClone());
  myclone->set_spos(std::move(cloneup), std::move(clonedn));
  return myclone;
}

} // namespace qmcplusplus
