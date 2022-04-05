//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
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

void SpinorSet::resetParameters(const opt_variables_type& optVariables){};

void SpinorSet::setOrbitalSetSize(int norbs) { OrbitalSetSize = norbs; };


void SpinorSet::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
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

void SpinorSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
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
                                 ValueVector& psi,
                                 GradVector& dpsi,
                                 ValueVector& d2psi,
                                 ValueVector& dspin)
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

void SpinorSet::mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                       const RefVectorWithLeader<ParticleSet>& P_list,
                                       int iat,
                                       const RefVector<ValueVector>& psi_v_list,
                                       const RefVector<GradVector>& dpsi_v_list,
                                       const RefVector<ValueVector>& d2psi_v_list,
                                       const RefVector<ValueVector>& dspin_v_list) const
{
  auto& spo_leader = spo_list.getCastedLeader<SpinorSet>();
  auto& P_leader   = P_list.getLeader();
  assert(this == &spo_leader);

  IndexType nw          = spo_list.size();
  SPOSet& up_spo_leader = *(spo_leader.spo_up);
  SPOSet& dn_spo_leader = *(spo_leader.spo_dn);
  RefVectorWithLeader<SPOSet> up_spo_list(up_spo_leader);
  RefVectorWithLeader<SPOSet> dn_spo_list(dn_spo_leader);
  up_spo_list.reserve(nw);
  dn_spo_list.reserve(nw);

  std::vector<ValueVector> mw_up_psi_work, mw_dn_psi_work;
  std::vector<GradVector> mw_up_dpsi_work, mw_dn_dpsi_work;
  std::vector<ValueVector> mw_up_d2psi_work, mw_dn_d2psi_work;
  mw_up_psi_work.reserve(nw);
  mw_up_dpsi_work.reserve(nw);
  mw_up_d2psi_work.reserve(nw);
  mw_dn_psi_work.reserve(nw);
  mw_dn_dpsi_work.reserve(nw);
  mw_dn_d2psi_work.reserve(nw);

  RefVector<ValueVector> up_psi_v_list, dn_psi_v_list;
  RefVector<GradVector> up_dpsi_v_list, dn_dpsi_v_list;
  RefVector<ValueVector> up_d2psi_v_list, dn_d2psi_v_list;
  up_psi_v_list.reserve(nw);
  up_dpsi_v_list.reserve(nw);
  up_d2psi_v_list.reserve(nw);
  dn_psi_v_list.reserve(nw);
  dn_dpsi_v_list.reserve(nw);
  dn_d2psi_v_list.reserve(nw);

  ValueVector tmp_val_vec(OrbitalSetSize);
  GradVector tmp_grad_vec(OrbitalSetSize);
  for (int iw = 0; iw < nw; iw++)
  {
    SpinorSet& spinor = spo_list.getCastedElement<SpinorSet>(iw);
    up_spo_list.emplace_back(*(spinor.spo_up));
    dn_spo_list.emplace_back(*(spinor.spo_dn));

    mw_up_psi_work.emplace_back(tmp_val_vec);
    up_psi_v_list.emplace_back(mw_up_psi_work.back());
    mw_dn_psi_work.emplace_back(tmp_val_vec);
    dn_psi_v_list.emplace_back(mw_dn_psi_work.back());

    mw_up_dpsi_work.emplace_back(tmp_grad_vec);
    up_dpsi_v_list.emplace_back(mw_up_dpsi_work.back());
    mw_dn_dpsi_work.emplace_back(tmp_grad_vec);
    dn_dpsi_v_list.emplace_back(mw_dn_dpsi_work.back());

    mw_up_d2psi_work.emplace_back(tmp_val_vec);
    up_d2psi_v_list.emplace_back(mw_up_d2psi_work.back());
    mw_dn_d2psi_work.emplace_back(tmp_val_vec);
    dn_d2psi_v_list.emplace_back(mw_dn_d2psi_work.back());
  }

  up_spo_leader.mw_evaluateVGL(up_spo_list, P_list, iat, up_psi_v_list, up_dpsi_v_list, up_d2psi_v_list);
  dn_spo_leader.mw_evaluateVGL(dn_spo_list, P_list, iat, dn_psi_v_list, dn_dpsi_v_list, dn_d2psi_v_list);

#pragma omp parallel for
  for (int iw = 0; iw < nw; iw++)
  {
    ParticleSet::Scalar_t s = P_list[iw].activeSpin(iat);
    RealType coss           = std::cos(s);
    RealType sins           = std::sin(s);

    ValueType eis(coss, sins);
    ValueType emis(coss, -sins);
    ValueType eye(0, 1.0);

    psi_v_list[iw].get()   = eis * up_psi_v_list[iw].get() + emis * dn_psi_v_list[iw].get();
    dpsi_v_list[iw].get()  = eis * up_dpsi_v_list[iw].get() + emis * dn_dpsi_v_list[iw].get();
    d2psi_v_list[iw].get() = eis * up_d2psi_v_list[iw].get() + emis * dn_d2psi_v_list[iw].get();
    dspin_v_list[iw].get() = eye * (eis * up_psi_v_list[iw].get() - emis * dn_psi_v_list[iw].get());
  }
}

void SpinorSet::evaluate_notranspose(const ParticleSet& P,
                                     int first,
                                     int last,
                                     ValueMatrix& logdet,
                                     GradMatrix& dlogdet,
                                     ValueMatrix& d2logdet)
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
                                        const RefVector<ValueMatrix>& logdet_list,
                                        const RefVector<GradMatrix>& dlogdet_list,
                                        const RefVector<ValueMatrix>& d2logdet_list) const
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

  std::vector<ValueMatrix> mw_up_logdet, mw_dn_logdet;
  std::vector<GradMatrix> mw_up_dlogdet, mw_dn_dlogdet;
  std::vector<ValueMatrix> mw_up_d2logdet, mw_dn_d2logdet;
  mw_up_logdet.reserve(nw);
  mw_dn_logdet.reserve(nw);
  mw_up_dlogdet.reserve(nw);
  mw_dn_dlogdet.reserve(nw);
  mw_up_d2logdet.reserve(nw);
  mw_dn_d2logdet.reserve(nw);

  RefVector<ValueMatrix> up_logdet_list, dn_logdet_list;
  RefVector<GradMatrix> up_dlogdet_list, dn_dlogdet_list;
  RefVector<ValueMatrix> up_d2logdet_list, dn_d2logdet_list;
  up_logdet_list.reserve(nw);
  dn_logdet_list.reserve(nw);
  up_dlogdet_list.reserve(nw);
  dn_dlogdet_list.reserve(nw);
  up_d2logdet_list.reserve(nw);
  dn_d2logdet_list.reserve(nw);

  ValueMatrix tmp_val_mat(nelec, OrbitalSetSize);
  GradMatrix tmp_grad_mat(nelec, OrbitalSetSize);
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

void SpinorSet::evaluate_notranspose_spin(const ParticleSet& P,
                                          int first,
                                          int last,
                                          ValueMatrix& logdet,
                                          GradMatrix& dlogdet,
                                          ValueMatrix& d2logdet,
                                          ValueMatrix& dspinlogdet)
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


void SpinorSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector& psi, ValueVector& dpsi)
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
