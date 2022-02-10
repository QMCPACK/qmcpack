//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-


#include "J2OrbitalSoA.h"
#include "CPU/SIMD/algorithm.hpp"
#include "BsplineFunctor.h"
#include "PadeFunctors.h"
#include "UserFunctor.h"
#include "FakeFunctor.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{
template<typename FT>
void J2OrbitalSoA<FT>::checkInVariables(opt_variables_type& active)
{
  for (auto& [key, functor] : J2Unique)
    functor->checkInVariables(active);
}

template<typename FT>
void J2OrbitalSoA<FT>::checkOutVariables(const opt_variables_type& active)
{
  myVars.clear();
  for (auto& [key, functor] : J2Unique)
  {
    functor->myVars.getIndex(active);
    myVars.insertFrom(functor->myVars);
  }
  // Remove inactive variables so the mappings are correct
  myVars.removeInactive();

  myVars.getIndex(active);

  const size_t NumVars = myVars.size();
  if (NumVars)
  {
    OffSet.resize(F.size());
    // Find first active variable for the starting offset
    int varoffset = -1;
    for (int i = 0; i < myVars.size(); i++)
    {
      varoffset = myVars.Index[i];
      if (varoffset != -1)
        break;
    }

    for (int i = 0; i < F.size(); ++i)
    {
      if (F[i] && F[i]->myVars.Index.size())
      {
        OffSet[i].first  = F[i]->myVars.Index.front() - varoffset;
        OffSet[i].second = F[i]->myVars.Index.size() + OffSet[i].first;
      }
      else
      {
        OffSet[i].first = OffSet[i].second = -1;
      }
    }
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::resetParameters(const opt_variables_type& active)
{
  if (!Optimizable)
    return;
  for (auto& [key, functor] : J2Unique)
    functor->resetParameters(active);
  for (int i = 0; i < myVars.size(); ++i)
  {
    int ii = myVars.Index[i];
    if (ii >= 0)
      myVars[i] = active[ii];
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::reportStatus(std::ostream& os)
{
  auto it(J2Unique.begin()), it_end(J2Unique.end());
  while (it != it_end)
  {
    (*it).second->myVars.print(os);
    ++it;
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  for (int k = 0; k < ratios.size(); ++k)
    ratios[k] =
        std::exp(Uat[VP.refPtcl] - computeU(VP.refPS, VP.refPtcl, VP.getDistTableAB(my_table_ID_).getDistRow(k)));
}

template<typename FT>
void J2OrbitalSoA<FT>::registerData(ParticleSet& P, WFBufferType& buf)
{
  if (Bytes_in_WFBuffer == 0)
  {
    Bytes_in_WFBuffer = buf.current();
    buf.add(Uat.begin(), Uat.end());
    buf.add(dUat.data(), dUat.end());
    buf.add(d2Uat.begin(), d2Uat.end());
    Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
    // free local space
    Uat.free();
    dUat.free();
    d2Uat.free();
  }
  else
  {
    buf.forward(Bytes_in_WFBuffer);
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  Uat.attachReference(buf.lendReference<valT>(N), N);
  dUat.attachReference(N, N_padded, buf.lendReference<valT>(N_padded * OHMMS_DIM));
  d2Uat.attachReference(buf.lendReference<valT>(N), N);
}

template<typename FT>
typename J2OrbitalSoA<FT>::LogValueType J2OrbitalSoA<FT>::updateBuffer(ParticleSet& P,
                                                                       WFBufferType& buf,
                                                                       bool fromscratch)
{
  evaluateGL(P, P.G, P.L, false);
  buf.forward(Bytes_in_WFBuffer);
  return log_value_;
}

template<typename FT>
typename J2OrbitalSoA<FT>::valT J2OrbitalSoA<FT>::computeU(const ParticleSet& P, int iat, const DistRow& dist)
{
  valT curUat(0);
  const int igt = P.GroupID[iat] * NumGroups;
  for (int jg = 0; jg < NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt + jg]);
    int iStart = P.first(jg);
    int iEnd   = P.last(jg);
    curUat += f2.evaluateV(iat, iStart, iEnd, dist.data(), DistCompressed.data());
  }
  return curUat;
}

template<typename FT>
typename J2OrbitalSoA<FT>::posT J2OrbitalSoA<FT>::accumulateG(const valT* restrict du, const DisplRow& displ) const
{
  posT grad;
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
  {
    const valT* restrict dX = displ.data(idim);
    valT s                  = valT();

#pragma omp simd reduction(+ : s) aligned(du, dX : QMC_SIMD_ALIGNMENT)
    for (int jat = 0; jat < N; ++jat)
      s += du[jat] * dX[jat];
    grad[idim] = s;
  }
  return grad;
}

template<typename FT>
J2OrbitalSoA<FT>::J2OrbitalSoA(const std::string& obj_name, ParticleSet& p)
    : WaveFunctionComponent("J2OrbitalSoA", obj_name),
      my_table_ID_(p.addTable(p, DTModes::NEED_TEMP_DATA_ON_HOST)),
      j2_ke_corr_helper(p, F)
{
  if (myName.empty())
    throw std::runtime_error("J2OrbitalSoA object name cannot be empty!");
  init(p);
  KEcorr = 0.0;
}

template<typename FT>
J2OrbitalSoA<FT>::~J2OrbitalSoA() = default;

template<typename FT>
void J2OrbitalSoA<FT>::init(ParticleSet& p)
{
  N         = p.getTotalNum();
  N_padded  = getAlignedSize<valT>(N);
  NumGroups = p.groups();

  Uat.resize(N);
  dUat.resize(N);
  d2Uat.resize(N);
  cur_u.resize(N);
  cur_du.resize(N);
  cur_d2u.resize(N);
  old_u.resize(N);
  old_du.resize(N);
  old_d2u.resize(N);
  F.resize(NumGroups * NumGroups, nullptr);
  DistCompressed.resize(N);
  DistIndice.resize(N);
}

template<typename FT>
void J2OrbitalSoA<FT>::addFunc(int ia, int ib, std::unique_ptr<FT> j)
{
  assert(ia < NumGroups);
  assert(ib < NumGroups);
  if (ia == ib)
  {
    if (ia == 0) //first time, assign everything
    {
      int ij = 0;
      for (int ig = 0; ig < NumGroups; ++ig)
        for (int jg = 0; jg < NumGroups; ++jg, ++ij)
          if (F[ij] == nullptr)
            F[ij] = j.get();
    }
    else
      F[ia * NumGroups + ib] = j.get();
  }
  else
  {
    // a very special case, 1 particle of each type (e.g. 1 up + 1 down)
    // uu/dd/etc. was prevented by the builder
    if (N == NumGroups)
      for (int ig = 0; ig < NumGroups; ++ig)
        F[ig * NumGroups + ig] = j.get();
    // generic case
    F[ia * NumGroups + ib] = j.get();
    F[ib * NumGroups + ia] = j.get();
  }
  std::stringstream aname;
  aname << ia << ib;
  J2Unique[aname.str()] = std::move(j);
}

template<typename FT>
std::unique_ptr<WaveFunctionComponent> J2OrbitalSoA<FT>::makeClone(ParticleSet& tqp) const
{
  auto j2copy = std::make_unique<J2OrbitalSoA<FT>>(myName, tqp);
  std::map<const FT*, FT*> fcmap;
  for (int ig = 0; ig < NumGroups; ++ig)
    for (int jg = ig; jg < NumGroups; ++jg)
    {
      int ij = ig * NumGroups + jg;
      if (F[ij] == 0)
        continue;
      typename std::map<const FT*, FT*>::iterator fit = fcmap.find(F[ij]);
      if (fit == fcmap.end())
      {
        auto fc      = std::make_unique<FT>(*F[ij]);
        fcmap[F[ij]] = fc.get();
        j2copy->addFunc(ig, jg, std::move(fc));
      }
    }
  j2copy->KEcorr      = KEcorr;
  j2copy->Optimizable = Optimizable;

  j2copy->myVars.clear();
  j2copy->myVars.insertFrom(myVars);
  j2copy->OffSet = OffSet;

  return j2copy;
}

/** intenal function to compute \f$\sum_j u(r_j), du/dr, d2u/dr2\f$
 * @param P particleset
 * @param iat particle index
 * @param dist starting distance
 * @param u starting value
 * @param du starting first deriv
 * @param d2u starting second deriv
 */
template<typename FT>
void J2OrbitalSoA<FT>::computeU3(const ParticleSet& P,
                                 int iat,
                                 const DistRow& dist,
                                 RealType* restrict u,
                                 RealType* restrict du,
                                 RealType* restrict d2u,
                                 bool triangle)
{
  const int jelmax = triangle ? iat : N;
  constexpr valT czero(0);
  std::fill_n(u, jelmax, czero);
  std::fill_n(du, jelmax, czero);
  std::fill_n(d2u, jelmax, czero);

  const int igt = P.GroupID[iat] * NumGroups;
  for (int jg = 0; jg < NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt + jg]);
    int iStart = P.first(jg);
    int iEnd   = std::min(jelmax, P.last(jg));
    f2.evaluateVGL(iat, iStart, iEnd, dist.data(), u, du, d2u, DistCompressed.data(), DistIndice.data());
  }
  //u[iat]=czero;
  //du[iat]=czero;
  //d2u[iat]=czero;
}

template<typename FT>
typename J2OrbitalSoA<FT>::PsiValueType J2OrbitalSoA<FT>::ratio(ParticleSet& P, int iat)
{
  //only ratio, ready to compute it again
  UpdateMode = ORB_PBYP_RATIO;
  cur_Uat    = computeU(P, iat, P.getDistTableAA(my_table_ID_).getTempDists());
  return std::exp(static_cast<PsiValueType>(Uat[iat] - cur_Uat));
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  const auto& d_table = P.getDistTableAA(my_table_ID_);
  const auto& dist    = d_table.getTempDists();

  for (int ig = 0; ig < NumGroups; ++ig)
  {
    const int igt = ig * NumGroups;
    valT sumU(0);
    for (int jg = 0; jg < NumGroups; ++jg)
    {
      const FuncType& f2(*F[igt + jg]);
      int iStart = P.first(jg);
      int iEnd   = P.last(jg);
      sumU += f2.evaluateV(-1, iStart, iEnd, dist.data(), DistCompressed.data());
    }

    for (int i = P.first(ig); i < P.last(ig); ++i)
    {
      // remove self-interaction
      const valT Uself = F[igt + ig]->evaluate(dist[i]);
      ratios[i]        = std::exp(Uat[i] + Uself - sumU);
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::GradType J2OrbitalSoA<FT>::evalGrad(ParticleSet& P, int iat)
{
  return GradType(dUat[iat]);
}

template<typename FT>
typename J2OrbitalSoA<FT>::PsiValueType J2OrbitalSoA<FT>::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  computeU3(P, iat, P.getDistTableAA(my_table_ID_).getTempDists(), cur_u.data(), cur_du.data(), cur_d2u.data());
  cur_Uat = simd::accumulate_n(cur_u.data(), N, valT());
  DiffVal = Uat[iat] - cur_Uat;
  grad_iat += accumulateG(cur_du.data(), P.getDistTableAA(my_table_ID_).getTempDispls());
  return std::exp(static_cast<PsiValueType>(DiffVal));
}

template<typename FT>
void J2OrbitalSoA<FT>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  // get the old u, du, d2u
  const auto& d_table = P.getDistTableAA(my_table_ID_);
  computeU3(P, iat, d_table.getOldDists(), old_u.data(), old_du.data(), old_d2u.data());
  if (UpdateMode == ORB_PBYP_RATIO)
  { //ratio-only during the move; need to compute derivatives
    const auto& dist = d_table.getTempDists();
    computeU3(P, iat, dist, cur_u.data(), cur_du.data(), cur_d2u.data());
  }

  valT cur_d2Uat(0);
  const auto& new_dr    = d_table.getTempDispls();
  const auto& old_dr    = d_table.getOldDispls();
  constexpr valT lapfac = OHMMS_DIM - RealType(1);
#pragma omp simd reduction(+ : cur_d2Uat)
  for (int jat = 0; jat < N; jat++)
  {
    const valT du   = cur_u[jat] - old_u[jat];
    const valT newl = cur_d2u[jat] + lapfac * cur_du[jat];
    const valT dl   = old_d2u[jat] + lapfac * old_du[jat] - newl;
    Uat[jat] += du;
    d2Uat[jat] += dl;
    cur_d2Uat -= newl;
  }
  posT cur_dUat;
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
  {
    const valT* restrict new_dX    = new_dr.data(idim);
    const valT* restrict old_dX    = old_dr.data(idim);
    const valT* restrict cur_du_pt = cur_du.data();
    const valT* restrict old_du_pt = old_du.data();
    valT* restrict save_g          = dUat.data(idim);
    valT cur_g                     = cur_dUat[idim];
#pragma omp simd reduction(+ : cur_g) aligned(old_dX, new_dX, save_g, cur_du_pt, old_du_pt : QMC_SIMD_ALIGNMENT)
    for (int jat = 0; jat < N; jat++)
    {
      const valT newg = cur_du_pt[jat] * new_dX[jat];
      const valT dg   = newg - old_du_pt[jat] * old_dX[jat];
      save_g[jat] -= dg;
      cur_g += newg;
    }
    cur_dUat[idim] = cur_g;
  }
  log_value_ += Uat[iat] - cur_Uat;
  Uat[iat]   = cur_Uat;
  dUat(iat)  = cur_dUat;
  d2Uat[iat] = cur_d2Uat;
}

template<typename FT>
void J2OrbitalSoA<FT>::recompute(const ParticleSet& P)
{
  const auto& d_table = P.getDistTableAA(my_table_ID_);
  for (int ig = 0; ig < NumGroups; ++ig)
  {
    for (int iat = P.first(ig), last = P.last(ig); iat < last; ++iat)
    {
      computeU3(P, iat, d_table.getDistRow(iat), cur_u.data(), cur_du.data(), cur_d2u.data(), true);
      Uat[iat] = simd::accumulate_n(cur_u.data(), iat, valT());
      posT grad;
      valT lap(0);
      const valT* restrict u   = cur_u.data();
      const valT* restrict du  = cur_du.data();
      const valT* restrict d2u = cur_d2u.data();
      const auto& displ        = d_table.getDisplRow(iat);
      constexpr valT lapfac    = OHMMS_DIM - RealType(1);
#pragma omp simd reduction(+ : lap) aligned(du, d2u : QMC_SIMD_ALIGNMENT)
      for (int jat = 0; jat < iat; ++jat)
        lap += d2u[jat] + lapfac * du[jat];
      for (int idim = 0; idim < OHMMS_DIM; ++idim)
      {
        const valT* restrict dX = displ.data(idim);
        valT s                  = valT();
#pragma omp simd reduction(+ : s) aligned(du, dX : QMC_SIMD_ALIGNMENT)
        for (int jat = 0; jat < iat; ++jat)
          s += du[jat] * dX[jat];
        grad[idim] = s;
      }
      dUat(iat)  = grad;
      d2Uat[iat] = -lap;
// add the contribution from the upper triangle
#pragma omp simd aligned(u, du, d2u : QMC_SIMD_ALIGNMENT)
      for (int jat = 0; jat < iat; jat++)
      {
        Uat[jat] += u[jat];
        d2Uat[jat] -= d2u[jat] + lapfac * du[jat];
      }
      for (int idim = 0; idim < OHMMS_DIM; ++idim)
      {
        valT* restrict save_g   = dUat.data(idim);
        const valT* restrict dX = displ.data(idim);
#pragma omp simd aligned(save_g, du, dX : QMC_SIMD_ALIGNMENT)
        for (int jat = 0; jat < iat; jat++)
          save_g[jat] -= du[jat] * dX[jat];
      }
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::LogValueType J2OrbitalSoA<FT>::evaluateLog(const ParticleSet& P,
                                                                      ParticleSet::ParticleGradient& G,
                                                                      ParticleSet::ParticleLaplacian& L)
{
  return evaluateGL(P, G, L, true);
}

template<typename FT>
WaveFunctionComponent::LogValueType J2OrbitalSoA<FT>::evaluateGL(const ParticleSet& P,
                                                                 ParticleSet::ParticleGradient& G,
                                                                 ParticleSet::ParticleLaplacian& L,
                                                                 bool fromscratch)
{
  if (fromscratch)
    recompute(P);
  log_value_ = valT(0);
  for (int iat = 0; iat < N; ++iat)
  {
    log_value_ += Uat[iat];
    G[iat] += dUat[iat];
    L[iat] += d2Uat[iat];
  }

  return log_value_ = -log_value_ * 0.5;
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi)
{
  log_value_ = 0.0;
  const auto& d_ee(P.getDistTableAA(my_table_ID_));
  valT dudr, d2udr2;

  Tensor<valT, DIM> ident;
  grad_grad_psi = 0.0;
  ident.diagonal(1.0);

  for (int i = 1; i < N; ++i)
  {
    const auto& dist  = d_ee.getDistRow(i);
    const auto& displ = d_ee.getDisplRow(i);
    auto ig           = P.GroupID[i];
    const int igt     = ig * NumGroups;
    for (int j = 0; j < i; ++j)
    {
      auto r    = dist[j];
      auto rinv = 1.0 / r;
      auto dr   = displ[j];
      auto jg   = P.GroupID[j];
      auto uij  = F[igt + jg]->evaluate(r, dudr, d2udr2);
      log_value_ -= uij;
      auto hess = rinv * rinv * outerProduct(dr, dr) * (d2udr2 - dudr * rinv) + ident * dudr * rinv;
      grad_grad_psi[i] -= hess;
      grad_grad_psi[j] -= hess;
    }
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateDerivatives(ParticleSet& P,
                                           const opt_variables_type& active,
                                           std::vector<ValueType>& dlogpsi,
                                           std::vector<ValueType>& dhpsioverpsi)
{
  if (myVars.size() == 0)
    return;

  evaluateDerivativesWF(P, active, dlogpsi);
  bool recalculate(false);
  std::vector<bool> rcsingles(myVars.size(), false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (active.recompute(kk))
      recalculate = true;
    rcsingles[k] = true;
  }
  if (recalculate)
  {
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (rcsingles[k])
      {
        dhpsioverpsi[kk] = -RealType(0.5) * ValueType(Sum(lapLogPsi[k])) - ValueType(Dot(P.G, gradLogPsi[k]));
      }
    }
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateDerivativesWF(ParticleSet& P,
                                             const opt_variables_type& active,
                                             std::vector<ValueType>& dlogpsi)
{
  if (myVars.size() == 0)
    return;

  resizeWFOptVectors();

  bool recalculate(false);
  std::vector<bool> rcsingles(myVars.size(), false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (active.recompute(kk))
      recalculate = true;
    rcsingles[k] = true;
  }
  if (recalculate)
  {
    ///precomputed recalculation switch
    std::vector<bool> RecalcSwitch(F.size(), false);
    for (int i = 0; i < F.size(); ++i)
    {
      if (OffSet[i].first < 0)
      {
        // nothing to optimize
        RecalcSwitch[i] = false;
      }
      else
      {
        bool recalcFunc(false);
        for (int rcs = OffSet[i].first; rcs < OffSet[i].second; rcs++)
          if (rcsingles[rcs] == true)
            recalcFunc = true;
        RecalcSwitch[i] = recalcFunc;
      }
    }
    dLogPsi              = 0.0;
    const size_t NumVars = myVars.size();
    for (int p = 0; p < NumVars; ++p)
    {
      gradLogPsi[p] = 0.0;
      lapLogPsi[p]  = 0.0;
    }
    std::vector<TinyVector<RealType, 3>> derivs(NumVars);
    const auto& d_table = P.getDistTableAA(my_table_ID_);
    constexpr RealType cone(1);
    constexpr RealType lapfac(OHMMS_DIM - cone);
    const size_t n  = d_table.sources();
    const size_t ng = P.groups();
    for (size_t i = 1; i < n; ++i)
    {
      const size_t ig   = P.GroupID[i] * ng;
      const auto& dist  = d_table.getDistRow(i);
      const auto& displ = d_table.getDisplRow(i);
      for (size_t j = 0; j < i; ++j)
      {
        const size_t ptype = ig + P.GroupID[j];
        if (RecalcSwitch[ptype])
        {
          std::fill(derivs.begin(), derivs.end(), 0.0);
          if (!F[ptype]->evaluateDerivatives(dist[j], derivs))
            continue;
          RealType rinv(cone / dist[j]);
          PosType dr(displ[j]);
          for (int p = OffSet[ptype].first, ip = 0; p < OffSet[ptype].second; ++p, ++ip)
          {
            RealType dudr(rinv * derivs[ip][1]);
            RealType lap(derivs[ip][2] + lapfac * dudr);
            //RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
            PosType gr(dudr * dr);
            dLogPsi[p] -= derivs[ip][0];
            gradLogPsi[p][i] += gr;
            gradLogPsi[p][j] -= gr;
            lapLogPsi[p][i] -= lap;
            lapLogPsi[p][j] -= lap;
          }
        }
      }
    }
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (rcsingles[k])
      {
        dlogpsi[kk] = dLogPsi[k];
      }
      //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(lapLogPsi[ip])-Dot(P.G,gradLogPsi[ip]));
    }
  }
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateDerivRatios(const VirtualParticleSet& VP,
                                           const opt_variables_type& optvars,
                                           std::vector<ValueType>& ratios,
                                           Matrix<ValueType>& dratios)
{
  evaluateRatios(VP, ratios);
  if (myVars.size() == 0)
    return;

  bool recalculate(false);
  std::vector<bool> rcsingles(myVars.size(), false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (optvars.recompute(kk))
      recalculate = true;
    rcsingles[k] = true;
  }

  if (recalculate)
  {
    ///precomputed recalculation switch
    std::vector<bool> RecalcSwitch(F.size(), false);
    for (int i = 0; i < F.size(); ++i)
    {
      if (OffSet[i].first < 0)
      {
        // nothing to optimize
        RecalcSwitch[i] = false;
      }
      else
      {
        bool recalcFunc(false);
        for (int rcs = OffSet[i].first; rcs < OffSet[i].second; rcs++)
          if (rcsingles[rcs] == true)
            recalcFunc = true;
        RecalcSwitch[i] = recalcFunc;
      }
    }
    const size_t NumVars = myVars.size();
    std::vector<RealType> derivs_ref(NumVars);
    std::vector<RealType> derivs(NumVars);
    const auto& d_table = VP.getDistTableAB(my_table_ID_);
    const size_t n      = d_table.sources();
    const size_t nt     = VP.getTotalNum();
    for (size_t i = 0; i < n; ++i)
    {
      if (i == VP.refPtcl)
        continue;
      const size_t ptype = VP.refPS.GroupID[i] * VP.refPS.groups() + VP.refPS.GroupID[VP.refPtcl];
      if (!RecalcSwitch[ptype])
        continue;
      const auto dist_ref = i < VP.refPtcl ? VP.refPS.getDistTableAA(my_table_ID_).getDistRow(VP.refPtcl)[i]
                                           : VP.refPS.getDistTableAA(my_table_ID_).getDistRow(i)[VP.refPtcl];
      //first calculate the old derivatives VP.refPtcl.
      std::fill(derivs_ref.begin(), derivs_ref.end(), 0.0);
      F[ptype]->evaluateDerivatives(dist_ref, derivs_ref);
      for (size_t j = 0; j < nt; ++j)
      {
        std::fill(derivs.begin(), derivs.end(), 0.0);
        F[ptype]->evaluateDerivatives(d_table.getDistRow(j)[i], derivs);
        for (int ip = 0, p = F[ptype]->myVars.Index.front(); ip < F[ptype]->myVars.Index.size(); ++ip, ++p)
          dratios[j][p] += derivs_ref[ip] - derivs[ip];
      }
    }
  }
}

template class J2OrbitalSoA<BsplineFunctor<QMCTraits::RealType>>;
template class J2OrbitalSoA<PadeFunctor<QMCTraits::RealType>>;
template class J2OrbitalSoA<UserFunctor<QMCTraits::RealType>>;
template class J2OrbitalSoA<FakeFunctor<QMCTraits::RealType>>;

} // namespace qmcplusplus
