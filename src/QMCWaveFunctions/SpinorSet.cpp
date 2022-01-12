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
    throw std::runtime_error("SpinorSet::set_spos(...):  up and down SPO components have different sizes.");

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
  assert(basis_up == basis_dn);
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

void SpinorSet::mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                                        const RefVectorWithLeader<ParticleSet>& P_list,
                                        int first,
                                        int last,
                                        const RefVector<ValueMatrix_t>& logdet_list,
                                        const RefVector<GradMatrix_t>& dlogdet_list,
                                        const RefVector<ValueMatrix_t>& d2logdet_list) const
{
  auto& spo_leader = spo_list.getCastedLeader<SpinorSet>();
  auto& P_leader   = P_list.getLeader();
  assert(this == &spo_leader);

  IndexType nw    = spo_list.size();
  IndexType nelec = P_leader.getTotalNum();

  SPOSet& up_spo_leader = *(spo_leader.spo_up);
  SPOSet& dn_spo_leader = *(spo_leader.spo_dn);
  RefVectorWithLeader<SPOSet> up_spo_list(up_spo_leader);
  RefVectorWithLeader<SPOSet> dn_spo_list(dn_spo_leader);
  up_spo_list.reserve(nw);
  dn_spo_list.reserve(nw);
  up_spo_list.clear();
  dn_spo_list.clear();

  std::vector<ValueMatrix_t> mw_up_logdet, mw_dn_logdet;
  std::vector<GradMatrix_t> mw_up_dlogdet, mw_dn_dlogdet;
  std::vector<ValueMatrix_t> mw_up_d2logdet, mw_dn_d2logdet;
  mw_up_logdet.reserve(nw);
  mw_dn_logdet.reserve(nw);
  mw_up_dlogdet.reserve(nw);
  mw_dn_dlogdet.reserve(nw);
  mw_up_d2logdet.reserve(nw);
  mw_dn_d2logdet.reserve(nw);
  mw_up_logdet.clear();
  mw_dn_logdet.clear();
  mw_up_dlogdet.clear();
  mw_dn_dlogdet.clear();
  mw_up_d2logdet.clear();
  mw_dn_d2logdet.clear();

  RefVector<ValueMatrix_t> up_logdet_list, dn_logdet_list;
  RefVector<GradMatrix_t> up_dlogdet_list, dn_dlogdet_list;
  RefVector<ValueMatrix_t> up_d2logdet_list, dn_d2logdet_list;
  up_logdet_list.reserve(nw);
  dn_logdet_list.reserve(nw);
  up_dlogdet_list.reserve(nw);
  dn_dlogdet_list.reserve(nw);
  up_d2logdet_list.reserve(nw);
  dn_d2logdet_list.reserve(nw);
  up_logdet_list.clear();
  dn_logdet_list.clear();
  up_dlogdet_list.clear();
  dn_dlogdet_list.clear();
  up_d2logdet_list.clear();
  dn_d2logdet_list.clear();

  ValueMatrix_t tmp_val_mat(nelec, OrbitalSetSize);
  GradMatrix_t tmp_grad_mat(nelec, OrbitalSetSize);
  for (int iw = 0; iw < nw; iw++)
  {
    SpinorSet& spinor = spo_list.getCastedElement<SpinorSet>(iw);
    up_spo_list.emplace_back(*(spinor.spo_up));
    dn_spo_list.emplace_back(*(spinor.spo_dn));

    mw_up_logdet.emplace_back(tmp_val_mat);
    up_logdet_list.emplace_back(mw_up_logdet.back());
    mw_dn_logdet.emplace_back(tmp_val_mat);
    dn_logdet_list.emplace_back(mw_dn_logdet.back());

    mw_up_dlogdet.emplace_back(tmp_grad_mat);
    up_dlogdet_list.emplace_back(mw_up_dlogdet.back());
    mw_dn_dlogdet.emplace_back(tmp_grad_mat);
    dn_dlogdet_list.emplace_back(mw_dn_dlogdet.back());

    mw_up_d2logdet.emplace_back(tmp_val_mat);
    up_d2logdet_list.emplace_back(mw_up_d2logdet.back());
    mw_dn_d2logdet.emplace_back(tmp_val_mat);
    dn_d2logdet_list.emplace_back(mw_dn_d2logdet.back());
  }

  up_spo_leader.mw_evaluate_notranspose(up_spo_list, P_list, first, last, up_logdet_list, up_dlogdet_list,
                                        up_d2logdet_list);
  dn_spo_leader.mw_evaluate_notranspose(dn_spo_list, P_list, first, last, dn_logdet_list, dn_dlogdet_list,
                                        dn_d2logdet_list);

#pragma omp parallel for
  for (int iw = 0; iw < nw; iw++)
  {
    for (int iat = 0; iat < nelec; iat++)
    {
      ParticleSet::Scalar_t s = P_list[iw].activeSpin(iat);
      RealType coss           = std::cos(s);
      RealType sins           = std::sin(s);
      ValueType eis(coss, sins);
      ValueType emis(coss, -sins);

      for (int no = 0; no < OrbitalSetSize; no++)
      {
        logdet_list[iw].get()(iat, no) =
            eis * up_logdet_list[iw].get()(iat, no) + emis * dn_logdet_list[iw].get()(iat, no);
        dlogdet_list[iw].get()(iat, no) =
            eis * up_dlogdet_list[iw].get()(iat, no) + emis * dn_dlogdet_list[iw].get()(iat, no);
        d2logdet_list[iw].get()(iat, no) =
            eis * up_d2logdet_list[iw].get()(iat, no) + emis * dn_d2logdet_list[iw].get()(iat, no);
      }
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

std::unique_ptr<SPOSet> SpinorSet::makeClone() const
{
  auto myclone = std::make_unique<SpinorSet>();
  std::unique_ptr<SPOSet> cloneup(spo_up->makeClone());
  std::unique_ptr<SPOSet> clonedn(spo_dn->makeClone());
  myclone->set_spos(std::move(cloneup), std::move(clonedn));
  return myclone;
}

} // namespace qmcplusplus
