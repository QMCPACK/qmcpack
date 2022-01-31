//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                  Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                  Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                  Shiv Upadhyay, shivnupadhyay@gmail.com, University of Pittsburgh
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ONEBODYSPINJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_ONEBODYSPINJASTROW_OPTIMIZED_SOA_H
#include "Configuration.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "Utilities/qmc_common.h"
#include "Utilities/IteratorUtility.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include <map>
#include <numeric>

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief Specialization for one-body Jastrow function using multiple functors
 */
template<class FT>
struct J1Spin : public WaveFunctionComponent
{
  ///alias FuncType
  using FuncType = FT;
  ///type of each component U, dU, d2U;
  using valT = typename FT::real_type;
  ///element position type
  using posT = TinyVector<valT, OHMMS_DIM>;
  ///use the same container
  using DistRow  = DistanceTable::DistRow;
  using DisplRow = DistanceTable::DisplRow;

  using GradDerivVec  = ParticleAttrib<QTFull::GradType>;
  using ValueDerivVec = ParticleAttrib<QTFull::ValueType>;

  ///table index
  const int myTableID;
  ///number of ions
  const int Nions;
  ///number of electrons
  const int Nelec;
  /* the number of ion groups if ions in 'Ions' particleset are grouped by species. 0 otherwise.
   * 0 Use slow code path. >= 1 use the code path with ion grouping
   */
  const int NumGroups;
  // number of electron groups
  const int NumTargetGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;

  ///variables handled by this orbital
  opt_variables_type myVars;

  valT curAt;
  valT curLap;
  posT curGrad;

  ///\f$Vat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Vat;
  aligned_vector<valT> U, dU, d2U, d3U;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  Vector<posT> Grad;
  Vector<valT> Lap;
  ///container for the unique Jastrow functors (size NumGroups x NumTargetGroups)
  std::vector<std::unique_ptr<FT>> J1UniqueFunctors;

  std::vector<std::pair<int, int>> OffSet;
  Vector<RealType> dLogPsi;
  std::vector<GradDerivVec> gradLogPsi;
  std::vector<ValueDerivVec> lapLogPsi;

  void resizeWFOptVectors()
  {
    dLogPsi.resize(myVars.size());
    gradLogPsi.resize(myVars.size(), GradDerivVec(Nelec));
    lapLogPsi.resize(myVars.size(), ValueDerivVec(Nelec));
  }

  J1Spin(const std::string& obj_name, const ParticleSet& ions, ParticleSet& els)
      : WaveFunctionComponent("J1Spin", obj_name),
        myTableID(els.addTable(ions)),
        Nions(ions.getTotalNum()),
        Nelec(els.getTotalNum()),
        NumGroups(determineNumGroups(ions)),
        NumTargetGroups(determineNumGroups(els)),
        Ions(ions)
  {
    if (myName.empty())
      throw std::runtime_error("J1Spin object name cannot be empty!");
    initialize(els);
  }

  J1Spin(const J1Spin& rhs) = delete;

  J1Spin(const J1Spin& rhs, ParticleSet& tqp)
      : WaveFunctionComponent("J1Spin", rhs.myName),
        myTableID(rhs.myTableID),
        Nions(rhs.Nions),
        Nelec(rhs.Nelec),
        NumGroups(rhs.NumGroups),
        NumTargetGroups(rhs.NumTargetGroups),
        Ions(rhs.Ions)
  {
    Optimizable = rhs.Optimizable;
    initialize(tqp);
    for (int i = 0; i < NumGroups; i++)
      for (int j = 0; j < NumTargetGroups; j++)
        if (rhs.J1UniqueFunctors[i * NumTargetGroups + j])
        {
          auto fc = std::make_unique<FT>(*rhs.J1UniqueFunctors[i * NumTargetGroups + j].get());
          addFunc(i, std::move(fc), j);
        }
    myVars = rhs.myVars;
    OffSet = rhs.OffSet;
  }

  /* determine NumGroups which controls the use of optimized code path using ion groups or not */
  static int determineNumGroups(const ParticleSet& ions)
  {
    const int num_species = ions.getSpeciesSet().getTotalNum();
    if (num_species == 1)
      return 1;
    else if (num_species > 1 && !ions.isGrouped())
      return 0;
    else
      return num_species;
  }

  /* initialize storage */
  void initialize(ParticleSet& els)
  {
    J1UniqueFunctors.resize(NumGroups * NumTargetGroups);
    Vat.resize(Nelec);
    Grad.resize(Nelec);
    Lap.resize(Nelec);

    U.resize(Nions);
    dU.resize(Nions);
    d2U.resize(Nions);
    d3U.resize(Nions);
    DistCompressed.resize(Nions);
    DistIndice.resize(Nions);
  }

  void addFunc(int source_type, std::unique_ptr<FT> afunc, int target_type = -1)
  {
    // if target type is not specified J1UniqueFunctors[i*NumTargetGroups + (0:NumTargetGroups)] is assigned
    // if target type is specified J1UniqueFunctors[i*NumTargetGroups + j] is assigned
    assert(target_type < NumTargetGroups);
    if (target_type == -1)
    {
      for (int i = 0; i < Nions; i++)
        for (int j = 0; j < NumTargetGroups; j++)
        {
          auto igroup = Ions.getGroupID(i);
          if (igroup == source_type && J1UniqueFunctors[igroup * NumTargetGroups + j] == nullptr)
            J1UniqueFunctors[igroup * NumTargetGroups + j] = std::move(afunc);
        }
    }
    else
    {
      for (int i = 0; i < Nions; i++)
        for (int j = 0; j < NumTargetGroups; j++)
        {
          auto igroup = Ions.getGroupID(i);
          if (Ions.getGroupID(i) == source_type && j == target_type &&
              J1UniqueFunctors[i * NumTargetGroups + j] == nullptr)
            J1UniqueFunctors[igroup * Nelec + j] = std::move(afunc);
        }
    }
  }

  void recompute(const ParticleSet& P) override
  {
    const auto& d_ie(P.getDistTableAB(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      computeU3(P, iat, d_ie.getDistRow(iat));
      Vat[iat] = simd::accumulate_n(U.data(), Nions, valT());
      Lap[iat] = accumulateGL(dU.data(), d2U.data(), d_ie.getDisplRow(iat), Grad[iat]);
    }
  }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    return evaluateGL(P, G, L, true);
  }

  void evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi) override
  {
    const auto& d_ie(P.getDistTableAB(myTableID));
    valT dudr, d2udr2;

    Tensor<valT, DIM> ident;
    grad_grad_psi = 0.0;
    ident.diagonal(1.0);

    for (int iel = 0; iel < Nelec; ++iel)
    {
      const auto& dist  = d_ie.getDistRow(iel);
      const auto& displ = d_ie.getDisplRow(iel);
      for (int iat = 0; iat < Nions; iat++)
      {
        auto gid   = Ions.getGroupID(iat) * NumTargetGroups + P.getGroupID(iel);
        auto* func = J1UniqueFunctors[gid].get();
        if (func != nullptr)
        {
          RealType r    = dist[iat];
          RealType rinv = 1.0 / r;
          PosType dr    = displ[iat];
          func->evaluate(r, dudr, d2udr2);
          grad_grad_psi[iel] -= rinv * rinv * outerProduct(dr, dr) * (d2udr2 - dudr * rinv) + ident * dudr * rinv;
        }
      }
    }
  }

  PsiValueType ratio(ParticleSet& P, int iat) override
  {
    UpdateMode = ORB_PBYP_RATIO;
    curAt      = computeU(P, iat, P.getDistTableAB(myTableID).getTempDists());
    return std::exp(static_cast<PsiValueType>(Vat[iat] - curAt));
  }

  inline void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override
  {
    for (int k = 0; k < ratios.size(); ++k)
      ratios[k] =
          std::exp(Vat[VP.refPtcl] - computeU(VP.refPS, VP.refPtcl, VP.getDistTableAB(myTableID).getDistRow(k)));
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
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

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, std::vector<ValueType>& dlogpsi) override
  {
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
      const size_t NumVars = myVars.size();
      for (int p = 0; p < NumVars; ++p)
      {
        gradLogPsi[p] = 0.0;
        lapLogPsi[p]  = 0.0;
      }
      dLogPsi = 0.0;

      const auto& d_table = P.getDistTableAB(myTableID);
      std::vector<TinyVector<RealType, 3>> derivs(NumVars);

      constexpr RealType cone(1);
      constexpr RealType lapfac(OHMMS_DIM - cone);
      const size_t ns = d_table.sources();
      const size_t nt = P.getTotalNum();

      aligned_vector<int> iadj(nt);
      aligned_vector<RealType> dist(nt);
      std::vector<PosType> displ(nt);

      for (size_t i = 0; i < ns; ++i)
      {
        for (size_t j = 0; j < nt; ++j)
        {
          const auto functor_idx = Ions.getGroupID(i) * NumTargetGroups + P.getGroupID(j);
          const int first(OffSet[functor_idx].first);
          const int last(OffSet[functor_idx].second);
          bool recalcFunc(false);
          for (int rcs = first; rcs < last; rcs++)
            if (rcsingles[rcs] == true)
              recalcFunc = true;
          if (recalcFunc)
          {
            auto* func = J1UniqueFunctors[functor_idx].get();
            if (func == nullptr)
              continue;
            std::fill(derivs.begin(), derivs.end(), 0);
            auto dist = P.getDistTableAB(myTableID).getDistRow(j)[i];
            if (!func->evaluateDerivatives(dist, derivs))
              continue;
            RealType rinv(cone / dist);
            const PosType& dr = P.getDistTableAB(myTableID).getDisplRow(j)[i];;
            for (int p = first, ip = 0; p < last; ++p, ++ip)
            {
              dLogPsi[p] -= derivs[ip][0];
              RealType dudr(rinv * derivs[ip][1]);
              gradLogPsi[p][j] += dudr * dr;
              lapLogPsi[p][j] -= derivs[ip][2] + lapfac * dudr;
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
          dlogpsi[kk] = ValueType(dLogPsi[k]);
        }
      }
    }
  }


  inline valT computeU(const ParticleSet& P, int iat, const DistRow& dist)
  {
    valT curVat(0);
    if (NumGroups > 0)
    {
      for (int jg = 0; jg < NumGroups; ++jg)
      {
        auto gid = jg * NumTargetGroups + P.getGroupID(iat);
        if (J1UniqueFunctors[gid])
          curVat +=
              J1UniqueFunctors[gid]->evaluateV(-1, Ions.first(jg), Ions.last(jg), dist.data(), DistCompressed.data());
      }
    }
    else
    {
      for (int c = 0; c < Nions; ++c)
      {
        auto gid = Ions.getGroupID(c) * NumTargetGroups + P.getGroupID(iat);
        if (J1UniqueFunctors[gid])
          curVat += J1UniqueFunctors[gid]->evaluate(dist[c]);
      }
    }
    return curVat;
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override
  {
    const auto& dist = P.getDistTableAB(myTableID).getTempDists();
    curAt            = valT(0);
    if (NumGroups > 0)
    {
      for (int ig = 0; ig < NumGroups; ++ig)
      {
        for (int jg = 0; jg < NumTargetGroups; ++jg)
        {
          auto gid = ig * NumTargetGroups + jg;
          if (J1UniqueFunctors[gid] != nullptr)
            curAt +=
                J1UniqueFunctors[gid]->evaluateV(-1, Ions.first(ig), Ions.last(ig), dist.data(), DistCompressed.data());
        }
      }
    }
    else
    {
      for (int ig = 0; ig < Nions; ++ig)
      {
        for (int jg = 0; jg < NumTargetGroups; ++jg)
        {
          auto gid = Ions.getGroupID(ig) * NumTargetGroups + jg;
          if (J1UniqueFunctors[gid] != nullptr)
            curAt += J1UniqueFunctors[gid]->evaluate(dist[ig]);
        }
      }
    }

    for (int i = 0; i < Nelec; ++i)
      ratios[i] = std::exp(Vat[i] - curAt);
  }

  inline LogValueType evaluateGL(const ParticleSet& P,
                                 ParticleSet::ParticleGradient& G,
                                 ParticleSet::ParticleLaplacian& L,
                                 bool fromscratch = false) override
  {
    if (fromscratch)
      recompute(P);

    for (size_t iat = 0; iat < Nelec; ++iat)
      G[iat] += Grad[iat];
    for (size_t iat = 0; iat < Nelec; ++iat)
      L[iat] -= Lap[iat];
    return log_value_ = -simd::accumulate_n(Vat.data(), Nelec, valT());
  }

  /** compute gradient and lap
   * @return lap
   */
  inline valT accumulateGL(const valT* restrict du, const valT* restrict d2u, const DisplRow& displ, posT& grad) const
  {
    valT lap(0);
    constexpr valT lapfac = OHMMS_DIM - RealType(1);
    //#pragma omp simd reduction(+:lap)
    for (int jat = 0; jat < Nions; ++jat)
      lap += d2u[jat] + lapfac * du[jat];
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
    {
      const valT* restrict dX = displ.data(idim);
      valT s                  = valT();
      //#pragma omp simd reduction(+:s)
      for (int jat = 0; jat < Nions; ++jat)
        s += du[jat] * dX[jat];
      grad[idim] = s;
    }
    return lap;
  }

  /** compute U, dU and d2U
   * @param P quantum particleset
   * @param iat the moving particle
   * @param dist starting address of the distances of the ions wrt the iat-th particle
   */
  inline void computeU3(const ParticleSet& P, int iat, const DistRow& dist)
  {
    if (NumGroups > 0)
    { //ions are grouped
      constexpr valT czero(0);
      std::fill_n(U.data(), Nions, czero);
      std::fill_n(dU.data(), Nions, czero);
      std::fill_n(d2U.data(), Nions, czero);

      for (int jg = 0; jg < NumGroups; ++jg)
      {
        auto gid = NumTargetGroups * jg + P.getGroupID(iat);
        if (J1UniqueFunctors[gid])
          J1UniqueFunctors[gid]->evaluateVGL(-1, Ions.first(jg), Ions.last(jg), dist.data(), U.data(), dU.data(),
                                             d2U.data(), DistCompressed.data(), DistIndice.data());
      }
    }
    else
    {
      for (int c = 0; c < Nions; ++c)
      {
        auto gid = Ions.getGroupID(c) * NumTargetGroups + P.getGroupID(iat);
        if (J1UniqueFunctors[gid])
        {
          U[c] = J1UniqueFunctors[gid]->evaluate(dist[c], dU[c], d2U[c]);
          dU[c] /= dist[c];
        }
      }
    }
  }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   */
  GradType evalGrad(ParticleSet& P, int iat) override { return GradType(Grad[iat]); }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   *
   * Using getTempDists(). curAt, curGrad and curLap are computed.
   */
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    UpdateMode = ORB_PBYP_PARTIAL;

    computeU3(P, iat, P.getDistTableAB(myTableID).getTempDists());
    curLap = accumulateGL(dU.data(), d2U.data(), P.getDistTableAB(myTableID).getTempDispls(), curGrad);
    curAt  = simd::accumulate_n(U.data(), Nions, valT());
    grad_iat += curGrad;
    return std::exp(static_cast<PsiValueType>(Vat[iat] - curAt));
  }

  /** Rejected move. Nothing to do */
  inline void restore(int iat) override {}

  /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    if (UpdateMode == ORB_PBYP_RATIO)
    {
      computeU3(P, iat, P.getDistTableAB(myTableID).getTempDists());
      curLap = accumulateGL(dU.data(), d2U.data(), P.getDistTableAB(myTableID).getTempDispls(), curGrad);
    }

    log_value_ += Vat[iat] - curAt;
    Vat[iat]  = curAt;
    Grad[iat] = curGrad;
    Lap[iat]  = curLap;
  }


  inline void registerData(ParticleSet& P, WFBufferType& buf) override
  {
    if (Bytes_in_WFBuffer == 0)
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Vat.begin(), Vat.end());
      buf.add(Grad.begin(), Grad.end());
      buf.add(Lap.begin(), Lap.end());
      Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
      // free local space
      Vat.free();
      Grad.free();
      Lap.free();
    }
    else
    {
      buf.forward(Bytes_in_WFBuffer);
    }
  }

  inline LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override
  {
    evaluateGL(P, P.G, P.L, false);
    buf.forward(Bytes_in_WFBuffer);
    return log_value_;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override
  {
    Vat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    Grad.attachReference(buf.lendReference<posT>(Nelec), Nelec);
    Lap.attachReference(buf.lendReference<valT>(Nelec), Nelec);
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override
  {
    auto cloned_J1Spin = std::make_unique<J1Spin<FT>>(*this, tqp);
    return cloned_J1Spin;
  }

  /**@{ WaveFunctionComponent virtual functions that are not essential for the development */
  void reportStatus(std::ostream& os) override
  {
    for (auto& J1UniqueFunctor : J1UniqueFunctors)
      if (J1UniqueFunctor != nullptr)
        J1UniqueFunctor->myVars.print(os);
  }

  void checkInVariables(opt_variables_type& active) override
  {
    myVars.clear();
    for (auto& J1UniqueFunctor : J1UniqueFunctors)
    {
      if (J1UniqueFunctor != nullptr)
      {
        J1UniqueFunctor->checkInVariables(active);
        J1UniqueFunctor->checkInVariables(myVars);
      }
    }
  }
  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.clear();
    for (auto& J1UniqueFunctor : J1UniqueFunctors)
    {
      if (J1UniqueFunctor)
      {
        J1UniqueFunctor->myVars.getIndex(active);
        myVars.insertFrom(J1UniqueFunctor->myVars);
      }
    }
    myVars.getIndex(active);
    const size_t NumVars = myVars.size();
    myVars.print(std::cout);
    if (NumVars)
    {
      OffSet.resize(J1UniqueFunctors.size());
      // Find first active variable for the starting offset
      int varoffset = -1;
      for (int i = 0; i < myVars.size(); i++)
      {
        varoffset = myVars.Index[i];
        if (varoffset != -1)
          break;
      }

      for (int i = 0; i < J1UniqueFunctors.size(); ++i)
      {
        if (J1UniqueFunctors[i] && J1UniqueFunctors[i]->myVars.Index.size())
        {
          OffSet[i].first  = J1UniqueFunctors[i]->myVars.Index.front() - varoffset;
          OffSet[i].second = J1UniqueFunctors[i]->myVars.Index.size() + OffSet[i].first;
        }
        else
        {
          OffSet[i].first = OffSet[i].second = -1;
        }
      }
    }
    Optimizable = myVars.is_optimizable();
    for (auto& J1UniqueFunctor : J1UniqueFunctors)
      if (J1UniqueFunctor != nullptr)
        J1UniqueFunctor->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active) override
  {
    if (!Optimizable)
      return;
    for (auto& J1UniqueFunctor : J1UniqueFunctors)
      if (J1UniqueFunctor != nullptr)
        J1UniqueFunctor->resetParameters(active);

    for (int i = 0; i < myVars.size(); ++i)
    {
      int ii = myVars.Index[i];
      if (ii >= 0)
        myVars[i] = active[ii];
    }
  }
  /**@} */

  inline GradType evalGradSource(ParticleSet& P, ParticleSet& source, int isrc) override
  {
    GradType g_return(0.0);
    const auto& d_ie(P.getDistTableAB(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      const auto& dist  = d_ie.getDistRow(iat);
      const auto& displ = d_ie.getDisplRow(iat);
      RealType r        = dist[isrc];
      RealType rinv     = 1.0 / r;
      PosType dr        = displ[isrc];
      auto gid          = source.getGroupID(isrc) * NumTargetGroups + P.getGroupID(iat);
      if (J1UniqueFunctors[gid] != nullptr)
      {
        U[isrc] = J1UniqueFunctors[gid]->evaluate(dist[isrc], dU[isrc], d2U[isrc], d3U[isrc]);
        g_return -= dU[isrc] * rinv * dr;
      }
    }
    return g_return;
  }

  inline GradType evalGradSource(ParticleSet& P,
                                 ParticleSet& source,
                                 int isrc,
                                 TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                                 TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override
  {
    GradType g_return(0.0);
    const auto& d_ie(P.getDistTableAB(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      const auto& dist  = d_ie.getDistRow(iat);
      const auto& displ = d_ie.getDisplRow(iat);
      RealType r        = dist[isrc];
      RealType rinv     = 1.0 / r;
      PosType dr        = displ[isrc];
      auto gid          = source.getGroupID(isrc) * NumTargetGroups + P.getGroupID(iat);
      if (J1UniqueFunctors[gid] != nullptr)
        U[isrc] = J1UniqueFunctors[gid]->evaluate(dist[isrc], dU[isrc], d2U[isrc], d3U[isrc]);
      else
        throw std::runtime_error("J1OrbitalSoa::evaluateGradSource:  J1UniqueFunctors[gid]==nullptr");

      g_return -= dU[isrc] * rinv * dr;

      //The following terms depend only on the radial component r.  Thus,
      //we compute them and mix with position vectors to acquire the full
      //cartesian vector objects.
      valT grad_component = (d2U[isrc] - dU[isrc] * rinv);
      valT lapl_component = d3U[isrc] + 2 * rinv * grad_component;

      for (int idim = 0; idim < OHMMS_DIM; idim++)
      {
        grad_grad[idim][iat] += dr[idim] * dr * rinv * rinv * grad_component;
        grad_grad[idim][iat][idim] += rinv * dU[isrc];

        lapl_grad[idim][iat] -= lapl_component * rinv * dr[idim];
      }
    }
    return g_return;
  }
};


} // namespace qmcplusplus
#endif
